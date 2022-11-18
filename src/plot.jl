using ArgParse, JLD2, Random, LinearAlgebra
using MultivariateStats

using Gadfly, Compose, DataFrames, StatsBase, Statistics
import Cairo, Fontconfig    
import ArgParse: parse_item
#using Revise

import Plots.scatter
import Plots.plot!
import Plots.savefig
import Plots


using UMAP: umap

function huecolors(n::Integer=100;alpha=0.8,saturation=1,brightness=1,precolors=[(0,1,0,alpha)],sufcolors=[])
    """
    "Flat Spike Trains to SpikeTimes and Trials, optionally sort trial based on `sv`"

    # Note code here is lifted from
    # https://github.com/Experica/NeuroAnalysis.jl/blob/master/src/Visualization/Visualization.jl#L41
    The neuroAnalysis package is not compatible with this packages build requirements so I hacked in their code as approriate. 
    Credit to the team at NeuroAnalysis.jl
    """
    hc = [(360(i-1)/n,saturation,brightness,alpha) for i in 1:n]
    prepend!(hc,precolors)
    append!(hc,sufcolors)
    hc
end
function flatspiketrains(xs,sv=[])
    """
    "Flat Spike Trains to SpikeTimes and Trials, optionally sort trial based on `sv`"

    # Note code here is lifted from
    # https://github.com/Experica/NeuroAnalysis.jl/blob/master/src/Visualization/Visualization.jl#L41
    The neuroAnalysis package is not compatible with this packages build requirements so I hacked in their code as approriate. 
    Credit to the team at NeuroAnalysis.jl
    """
    tn = length(xs)
    if isempty(sv)
        issort=false
    elseif length(sv)==tn
        issort=true
    else
        @warn """Length of "xs" and "sv" do not match, sorting ignored."""
        issort=false
    end
    if issort
        sxs=xs[sortperm(sv)]
        ssv=sort(sv)
    else
        sxs=xs
        ssv=sv
    end
    x=Float64[];y=Float64[];s=[]
    for i in 1:tn
        v = sxs[i];n=length(v)
        n==0 && continue
        append!(x,v);append!(y,fill(i,n))
        issort && append!(s,fill(ssv[i],n))
    end
    return x,y,s
end

function vstack(xs)
    """
    "Flat Spike Trains to SpikeTimes and Trials, optionally sort trial based on `sv`"

    # Note code here is lifted from
    # https://github.com/Experica/NeuroAnalysis.jl/blob/master/src/Visualization/Visualization.jl#L41
    The neuroAnalysis package is not compatible with this packages build requirements so I hacked in their code as approriate. 
    Credit to the team at NeuroAnalysis.jl
    """
    tn = length(xs)
    n = length(xs[1])
    mat = Matrix{Float64}(undef,tn,n)
    for i in 1:tn
        mat[i,:] = xs[i]
    end
    return mat
end

#using Plots,StatsPlots,VegaLite
#import Plots: cgrad
"scatter plot of spike trains"
function plotspiketrain(x,y;group::Vector=[],timeline=[0],color=huecolors(length(unique(group))),title="",size=(800,600))
    """
    "Flat Spike Trains to SpikeTimes and Trials, optionally sort trial based on `sv`"

    # Note code here is lifted from
    # https://github.com/Experica/NeuroAnalysis.jl/blob/master/src/Visualization/Visualization.jl#L41
    The neuroAnalysis package is not compatible with this packages build requirements so I hacked in their code as approriate. 
    Credit to the team at NeuroAnalysis.jl
    """
    nt = isempty(x) ? 0 : maximum(y)
    ms =0.45*size[2]/(nt+5)
    p = Plots.plot(;size,leg=false,title,grid=false)
    if isempty(group)
        Plots.scatter!(p,x,y;label="SpikeTrain",markershape=:vline,markersize=ms,markerstrokewidth = 0)
    else
        Plots.scatter!(p,x,y;group,markershape=:vline,markersize=ms,markerstrokewidth = 0,color=permutedims(color))
    end
    Plots.vline!(p,timeline;line=(:grey),label="TimeLine",xaxis="ms",yaxis=("Trial"))
end
function plotspiketrain(sts::Vector;uids::Vector=[],sortvalues=[],timeline=[0],color=huecolors(0),title="",size=(800,600))
    if isempty(uids)
        g=uids;uc=color
    else
        fuids = flatspiketrains(uids,sortvalues)[1]
        #@show(fuids)
        g=map(i->"U$i",fuids);uc=huecolors(length(unique(fuids)))
    end
    plotspiketrain(flatspiketrains(sts,sortvalues)[1:2]...;group=g,timeline,color=uc,title,size)
end


function ArgParse.parse_item(::Type{Vector{Int}}, x::AbstractString)
    return eval(Meta.parse(x))

end

if !(@isdefined nss)
    s = ArgParseSettings()
    @add_arg_table! s begin
    "--ineurons_to_plot", "-i"
        help = "which neurons to plot.  must be the same or a subset of ineurons_to_test used in test.jl"
        arg_type = Vector{Int}
    "test_file"
        help = "full path to the JLD file output by test.jl.  this same directory needs to contain the parameters in param.jld2, the synaptic targets in xtarg.jld2, and (optionally) the spike rate in rate.jld2"
        required = true
    end

    parsed_args = parse_args(s)
    #@show(parsed_args["test_file"])
    d = load(parsed_args["test_file"])
    ineurons_to_test = d["ineurons_to_test"]
    @show(length(ineurons_to_test))
    ineurons_to_plot = something(parsed_args["ineurons_to_plot"], ineurons_to_test)
    @show(length(ineurons_to_plot))

    all(in.(ineurons_to_plot,[ineurons_to_test])) || error("ineurons_to_plot must be the same or a subset of ineurons_to_test in test.jl")

    if all(in.(ineurons_to_test,[ineurons_to_plot]))
        nss = d["nss"]
        #@show(nss)

        timess = d["timess"]
        xtotals = d["xtotals"]
    else
        itest = []
        for iplot in ineurons_to_plot
            push!(itest, findfirst(iplot .== ineurons_to_test))
        end
        nss = similar(d["nss"])
        timess = similar(d["timess"])
        xtotals = similar(d["xtotals"])
        for ij in eachindex(d["nss"])
            nss[ij] = d["nss"][ij][itest]
            timess[ij] = d["timess"][ij][itest,:]
            xtotals[ij] = d["xtotals"][ij][:,itest]
        end
    end
    tt = convert(Array{Float64,2},timess[1])
    println("Number of spikes")
    println(length(tt[tt .> 0]))
    whole_mat = reduce(vcat, nss)
    plotspiketrain(whole_mat)#|>display
    savefig("Spike_Raster2.png")

    M = fit(PCA, tt; maxoutdim=2)
    Yte = predict(M, tt)


    res_jl = umap(tt; n_neighbors=5, min_dist=0.001, n_epochs=100)
    println(length(res_jl[1,:]))
    scatter(res_jl[1,:], res_jl[2,:], 
            title="spike train: Julia UMAP", marker=(2, 2, :auto))#|>display
    savefig("UMAP2.png")

    #p = scatter(Yte[1,:],Yte[2,:],Yte[3,:],marker=:circle,linewidth=1)
    #plot!(p,xlabel="PC1",ylabel="PC2",zlabel="PC3")|>display

    p = scatter(Yte[1,:],Yte[2,:],marker=:circle,linewidth=1)
    plot!(p,xlabel="PC1",ylabel="PC2")#|>display

    plot!(p)
    savefig("PCA2.png")



    include(joinpath(@__DIR__,"struct.jl"))
    p = load(joinpath(dirname(parsed_args["test_file"]),"param.jld2"), "p")

    xtarg = load(joinpath(dirname(parsed_args["test_file"]),"xtarg.jld2"), "xtarg")
    if isfile(joinpath(dirname(parsed_args["test_file"]),"rate.jld2"))
        rate = load(joinpath(dirname(parsed_args["test_file"]),"rate.jld2"), "rate")
    else
        rate = missing
    end

    output_prefix = splitext(parsed_args["test_file"])[1]
else
    ineurons_to_plot = parsed_args["ineurons_to_test"]

    xtarg = load(joinpath(parsed_args["data_dir"],"xtarg.jld2"), "xtarg")
    if isfile(joinpath(parsed_args["data_dir"],"rate.jld2"))
        rate = load(joinpath(parsed_args["data_dir"],"rate.jld2"), "rate")
    else
        rate = missing
    end

    output_prefix = joinpath(parsed_args["data_dir"], "test")
end


ntrials = size(nss,1)
ntasks = size(nss,2)
nneurons = length(nss[1])
nrows = isqrt(nneurons)
ncols = cld(nneurons, nrows)
@show(ncols)
@show(nrows)

for itask = 1:ntasks

ps = Union{Plot,Context}[]
for ci=1:nneurons
    df = DataFrame((t = (1:size(xtarg,1)).*p.learn_every/1000,
                    xtarg = xtarg[:,ineurons_to_plot[ci],itask],
                    xtotal1 = xtotals[1,itask][:,ci]))
    xtotal_ci = hcat((x[:,ci] for x in xtotals[:,itask])...)
    df[!,:xtotal_mean] = dropdims(mean(xtotal_ci, dims=2), dims=2)
    df[!,:xtotal_std] = dropdims(std(xtotal_ci, dims=2), dims=2)
    transform!(df, [:xtotal_mean, :xtotal_std] => ByRow((mu,sigma)->mu+sigma) => :xtotal_upper)
    transform!(df, [:xtotal_mean, :xtotal_std] => ByRow((mu,sigma)->mu-sigma) => :xtotal_lower)
    push!(ps, plot(df, x=:t, y=Col.value(:xtarg, :xtotal_mean, :xtotal1),
                   color=Col.index(:xtarg, :xtotal_mean, :xtotal1),
                   ymax=Col.value(:xtotal_upper), ymin=Col.value(:xtotal_lower),
                   Geom.line, Geom.ribbon,
                   Guide.colorkey(title="", labels=["data","model","model1"]),
                   Guide.title(string("neuron #", ineurons_to_plot[ci])),
                   Guide.xlabel("time (sec)", orientation=:horizontal),
                   Guide.ylabel("synaptic input", orientation=:vertical),
                   Guide.xticks(orientation=:horizontal)))
end
append!(ps, fill(Compose.context(), nrows*ncols-nneurons))
gridstack(permutedims(reshape(ps, ncols, nrows), (2,1))) |>
        PDF(string(output_prefix, "-syninput-task$itask.pdf"), 8cm*ncols, 6.5cm*nrows)

timess_cat = hcat(timess[:,itask]...)
ps = Union{Plot,Context}[]
for ci=1:nneurons
    psth = fit(Histogram, vec(timess_cat[ci,:]), p.stim_off : p.learn_every : p.train_time)
    df = DataFrame(t=p.learn_every/1000 : p.learn_every/1000 : p.train_time/1000-1,
                   model=psth.weights./ntrials./p.learn_every*1000)
    if ismissing(rate)
        scale_color = Scale.color_discrete(n->Scale.default_discrete_colors(n+1)[2:end])
        cols = (:model, )
    else
        scale_color = Scale.default_discrete_colors
        df[!,:data] = rate[:, ineurons_to_plot[ci]]
        cols = (:data, :model)
    end
    push!(ps, plot(df, x=:t, y=Col.value(cols...), color=Col.index(cols...),
                   Geom.line,
                   scale_color,
                   Guide.colorkey(title=""),
                   Guide.title(string("neuron #", ineurons_to_plot[ci])),
                   Guide.xlabel("time (sec)", orientation=:horizontal),
                   Guide.ylabel("spike rate", orientation=:vertical),
                   Guide.xticks(orientation=:horizontal)))
end
append!(ps, fill(Compose.context(), nrows*ncols-nneurons))
gridstack(permutedims(reshape(ps, ncols, nrows), (2,1))) |>
        PDF(string(output_prefix , "-psth-task$itask.pdf"), 8cm*ncols, 6.5cm*nrows)

end
