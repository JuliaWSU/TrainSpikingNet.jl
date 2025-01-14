using LinearAlgebra, Random, JLD2, Statistics, StatsBase, ArgParse, SymmetricFormats
using Plots
using StatsBase
using Random
using ProgressMeter
using StatsPlots
using Plots
#using UnicodePlots
unicodeplots()
#using Debugger
using JuliaInterpreter
function ArgParse.parse_item(::Type{Vector{Int}}, x::AbstractString)
    return eval(Meta.parse(x))
end
#verbose = false
verbose = true

aps = ArgParseSettings()

@add_arg_table! aps begin
    "data_dir"
        help = "full path to the directory containing the parameters file"
        required = true
    "--itasks", "-t"
        help = "a vector specifying which tasks to learn"
        arg_type = Vector{Int}
        default = [1]
        range_tester = x->all(x.>0)
end

add_arg_group!(aps, "mutually exclusive arguments.  if neither is specified, sinusoids\nwill be generated for synpatic inputs", exclusive = true);

@add_arg_table! aps begin
    "--xtarg_file", "-x"
        help = "full path to the JLD file containing the synaptic current targets"
    "--spikerate_file", "-s"
        help = "full path to the JLD file containing the spike rates"
end

parsed_args = parse_args(aps)

# --- set up variables --- #
include(joinpath(@__DIR__,"struct.jl"))
include(joinpath(parsed_args["data_dir"],"param.jl"))
include(joinpath(@__DIR__,"cpu","variables.jl"))

if p.Ncells == typemax(IntPrecision)
  @warn "IntPrecision is too small for GPU (but fine for CPU)"
elseif p.Ncells > typemax(IntPrecision)
  @error "IntPrecision is too small"
end

# --- load code --- #
macro maybethread(loop)
  if Threads.nthreads()>1
    quote Threads.@threads $(Expr(loop.head,
                             Expr(loop.args[1].head, esc.(loop.args[1].args)...),
                             esc(loop.args[2]))); end
  else
    @warn "running single threaded"
    quote $(esc(loop)); end 
  end
end

kind=:init
include(p.genStaticWeights_file)
include(p.genPlasticWeights_file)
include(p.genFfwdRate_file)
include(p.genTarget_file)
include(p.genStim_file)
include(joinpath("cpu","loop.jl"))
include("rate2synInput.jl")

#----------- initialization --------------#
w0Index, w0Weights, nc0 = genStaticWeights(p.genStaticWeights_args)
UnicodePlots.spy(w0Weights) |> display
#UnicodePlots.spy(w0Index) |> display

ffwdRate = genFfwdRate(p.genFfwdRate_args)

itask = 1

times = []
for i in 1:p.Ncells
    push!(times,Float32[])
end
uavg, ns0, ustd,times = loop_init(itask, nothing, nothing, p.stim_off, p.train_time, dt,
    p.Nsteps, p.Ncells, p.Ne, nothing, refrac, vre, invtauedecay,
    invtauidecay, nothing, mu, thresh, tau, nothing, times, ns, nothing,
    ns_ffwd, forwardInputsE, forwardInputsI, nothing, forwardInputsEPrev,
    forwardInputsIPrev, nothing, nothing, nothing, xedecay, xidecay, nothing,
    synInputBalanced, synInput, nothing, nothing, bias, nothing, nothing, lastSpike,
    nothing, nothing, nothing, nothing, nothing, v, rng, noise, rndFfwd, sig,
    nothing, nothing, w0Index, w0Weights, nc0, nothing, nothing, nothing,
    nothing, nothing, nothing, nothing, nothing, nothing, nothing, uavg,
    utmp, ffwdRate)



function convert_ragged_arraytodense(times)
    converttimes = []
    convertnodes = []
    for (ind,i) in enumerate(times)
        for t in i
            append!(converttimes,i)
            append!(convertnodes,ind)

        end

    end
    (converttimes,convertnodes)
end
    

(times,nodes) = convert_ragged_arraytodense(times)

wpWeightFfwd, wpWeightIn, wpIndexIn, ncpIn =
    genPlasticWeights(p.genPlasticWeights_args, w0Index, nc0, ns0)
if verbose
    UnicodePlots.spy(wpWeightIn) |> display
    #UnicodePlots.spy(wpIndexIn) |> display

end


# get indices of postsynaptic cells for each presynaptic cell
wpIndexConvert = zeros(Int, p.Ncells,p.Ncells)# p.Lexc+p.Linh)
wpIndexOutD = Dict()
ncpOut = Array{Int}(undef, p.Ncells)
for i = 1:p.Ncells
    wpIndexOutD[i] = []
end


for postCell = 1:p.Ncells
    if postCell!=0
        for i = 1:ncpIn[postCell]
 
            ##
            # wpIndexIn
            # incompatible with indexs provided.
            ##
            if postCell <= sizeof(wpIndexIn[:,i]) # unfortunate hack
                # to do figure out how to get rid of unfortunate hack.
                preCell = wpIndexIn[postCell,i] # This line fails
                if preCell!=0

                    push!(wpIndexOutD[preCell], postCell)

                    wpIndexConvert[postCell,i] = length(wpIndexOutD[preCell])
                end
            end
        end
    end
end
for preCell = 1:p.Ncells
    if preCell!=0


        ncpOut[preCell] = length(wpIndexOutD[preCell])
    end
end



# get weight, index of outgoing connections
wpIndexOut = zeros(Int, maximum(ncpOut),p.Ncells)
for preCell = 1:p.Ncells-1

    wpIndexOut[1:ncpOut[preCell],preCell] = wpIndexOutD[preCell]
end


Ntime = floor(Int, (p.train_time-p.stim_off)/p.learn_every)
if parsed_args["xtarg_file"] !== nothing
    xtarg_dict = load(parsed_args["xtarg_file"])
    xtarg = xtarg_dict[first(keys(xtarg_dict))]
    if size(xtarg,1) != Ntime
        error(parsed_args["xtarg_file"],
              " should have (train_time-stim_off)/learn_every = ",
              Ntime, " rows")
    end
    ndims(xtarg)==2 && (xtarg = xtarg[:,:,[CartesianIndex()]])
    if any(parsed_args["itasks"] .> size(xtarg,3))
        error("an element of --itasks exceeds the size of the third dimension of ",
              parsed_args["xtarg_file"])
    end
    xtarg = xtarg[:,:,parsed_args["itasks"]]
elseif parsed_args["spikerate_file"] !== nothing
    xtarg = rate2synInput(p, p.K==0 ? sig : (ustd / sqrt(p.tauedecay * 1.3))) # factor 1.3 was calibrated manually
else
    xtarg = Array{Float64}(undef, Ntime, p.Ncells, length(parsed_args["itasks"]))
    for itask = 1:length(parsed_args["itasks"])
        xtarg[:,:,itask] = genTarget(p.genTarget_args, uavg)
    end
end
timeSteps = round(Int, (p.stim_off - p.stim_on) / dt)
stim = Array{Float64}(undef, timeSteps, p.Ncells, length(parsed_args["itasks"]))
for itask = 1:length(parsed_args["itasks"])
    stim[:,:,itask] = genStim(p.genStim_args)
end

# --- set up correlation matrix --- #
pLrec = p.Lexc + p.Linh;

# L2-penalty
Pinv_L2 = Diagonal(repeat([p.penlambda], pLrec));
# row sum penalty
vec10 = [ones(p.Lexc); zeros(p.Linh)];
vec01 = [zeros(p.Lexc); ones(p.Linh)];
Pinv_rowsum = p.penmu*(vec10*vec10' + vec01*vec01');
# sum of penalties
Pinv_rec = Pinv_L2 + Pinv_rowsum;
Pinv_ffwd = Diagonal(repeat([p.penlamFF], p.Lffwd))
Pinv = zeros(pLrec+p.Lffwd, pLrec+p.Lffwd)
Pinv[1:pLrec, 1:pLrec] = Pinv_rec
Pinv[pLrec+1 : pLrec+p.Lffwd, pLrec+1 : pLrec+p.Lffwd] = Pinv_ffwd
Pinv_norm = p.PType(Symmetric(UpperTriangular(Pinv) \ I));


println("got here")
#----------- save initialization --------------#
save(joinpath(parsed_args["data_dir"],"param.jld2"), "p", p)
save(joinpath(parsed_args["data_dir"],"w0Index.jld2"), "w0Index", w0Index)
save(joinpath(parsed_args["data_dir"],"w0Weights.jld2"), "w0Weights", w0Weights)
save(joinpath(parsed_args["data_dir"],"nc0.jld2"), "nc0", nc0)
save(joinpath(parsed_args["data_dir"],"stim.jld2"), "stim", stim)
save(joinpath(parsed_args["data_dir"],"xtarg.jld2"), "xtarg", xtarg)
save(joinpath(parsed_args["data_dir"],"wpIndexIn.jld2"), "wpIndexIn", wpIndexIn)
save(joinpath(parsed_args["data_dir"],"wpIndexOut.jld2"), "wpIndexOut", wpIndexOut)
save(joinpath(parsed_args["data_dir"],"wpIndexConvert.jld2"), "wpIndexConvert", wpIndexConvert)
save(joinpath(parsed_args["data_dir"],"wpWeightFfwd.jld2"), "wpWeightFfwd", wpWeightFfwd)
save(joinpath(parsed_args["data_dir"],"wpWeightIn.jld2"), "wpWeightIn", wpWeightIn)
save(joinpath(parsed_args["data_dir"],"ncpIn.jld2"), "ncpIn", ncpIn)
save(joinpath(parsed_args["data_dir"],"ncpOut.jld2"), "ncpOut", ncpOut)
save(joinpath(parsed_args["data_dir"],"P.jld2"), "P", Pinv_norm)
save(joinpath(parsed_args["data_dir"],"ffwdRate.jld2"), "ffwdRate", ffwdRate)
