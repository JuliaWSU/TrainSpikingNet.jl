PPrecision = Float32  # can be <:Integer on GPUs
PScale = 1  # if PPrecision<:Integer then PScale should be e.g. 2^(nbits-2)
FloatPrecision = Float32
IntPrecision = UInt16
PType=Symmetric  # or SymmetricPacked for large models
seed=1

example_neurons = 25  # no. of neurons to save for visualization 
wid = 50  # width (ms) of the moving average window in time

monitor_resources_used = 0  # set to N to measure every N seconds
 
rng_func = Dict("gpu"=>:(MersenneTwister()), "cpu"=>:(MersenneTwister()))
rng = eval(rng_func["cpu"])
isnothing(seed) || Random.seed!(rng, seed)
save(joinpath(parsed_args["data_dir"],"rng-init.jld2"), "rng", rng)

dt = 0.1 #simulation timestep (ms)

# network size
Ncells = 1024
Ne = floor(Int, Ncells*0.5)
Ni = ceil(Int, Ncells*0.5)

# innate, train, test time (ms)
train_duration = 1000.0
stim_on        = 800.0
stim_off       = 1000.0
train_time     = stim_off + train_duration

Nsteps = round(Int, train_time/dt)

choose_task_func = :((iloop, ntasks) -> iloop % ntasks + 1)   # or rand(1:ntasks)

genStim_file = "genStim-random.jl"
genStim_args = Dict(:stim_on => stim_on, :stim_off => stim_off, :dt => dt, :Ncells => Ncells,
                    :mu => 0.0, :b => 1/20, :sig => 0.2)

# training variables
penlambda      = 3.0 # 0.1 or 0.5 
penlamFF       = 1.0
penmu          = 8.0 # 2.0
frac           = 1.0
learn_every    = 10.0 # (ms)

genTarget_file = "genTarget-sinusoids.jl"
genTarget_args = Dict(:train_time => train_time, :stim_off => stim_off, :learn_every => learn_every, :Ncells => Ncells, :Nsteps => Nsteps, :dt => dt,
                      :A => 0.5, :period => 1000.0, :biasType => "zero", :mu_ou_bias => 0.0, :b_ou_bias => 1/400, :sig_ou_bias => 0.02)

# neuron param      
taue = 10 #membrane time constant for exc. neurons (ms)
taui = 10 
threshe = 1.0 # spike threshold
threshi = 1.0   
refrac = 0.1 # refractory period
vre = 0.0

genStaticWeights_file = "genStaticWeights-erdos-renyi.jl"
genStaticWeights_args = Dict(:Ncells => Ncells, :Ne => Ne,
                             :pree => 0.1, :prie => 0.1, :prei => 0.1, :prii => 0.1)

K = round(Int, Ne*genStaticWeights_args[:pree])
sqrtK = sqrt(K)

g = 1.0
je = 2.0 / sqrtK * taue * g
ji = 2.0 / sqrtK * taue * g 
jx = 0.08 * sqrtK * g 

merge!(genStaticWeights_args, Dict(:jee => 0.15je, :jie => je, :jei => -0.75ji, :jii => -ji))

# plastic weights
L = round(Int,sqrt(K)*2.0) # number of exc/inh plastic weights per neuron
Lffwd = 0
Lexc = L # excitatory L
Linh = L # inhibitory L

genFfwdRate_file = "genFfwdRate-random.jl"
genFfwdRate_args = Dict(:train_time => train_time, :stim_off => stim_off, :dt => dt, :Lffwd => Lffwd,
                        :mu => 5, :bou => 1/400, :sig => 0.2, :wid => 500)

#synaptic time constants (ms) 
tauedecay = 3
tauidecay = 3
taudecay_plastic = 150  # can be a vector too, e.g. (150-70)*rand(rng, Ncells) .+ 70

muemin = jx*1.5 # exc external input
muemax = jx*1.5
muimin = jx # inh external input
muimax = jx

mu = Vector{Float64}(undef, Ncells)
mu[1:Ne] = (muemax-muemin)*rand(rng, Ne) .+ muemin
mu[(Ne+1):Ncells] = (muimax-muimin)*rand(rng, Ni) .+ muimin

wpscale = sqrt(L) * 2.0

genPlasticWeights_file = "genPlasticWeights-erdos-renyi.jl"
genPlasticWeights_args = Dict(:Ncells => Ncells, :frac => frac, :Ne => Ne, :L => L, :Lexc => Lexc, :Linh => Linh, :Lffwd => Lffwd,
                              :wpee => 2.0 * taue * g / wpscale,
                              :wpie => 2.0 * taue * g / wpscale,
                              :wpei => -2.0 * taue * g / wpscale,
                              :wpii => -2.0 * taue * g / wpscale,
                              :wpffwd => 0)

noise_model=:current  # or :voltage
sig = 0  # 0.65

correlation_var = K>0 ? :xtotal : :xplastic

maxrate = 500 #(Hz) maximum average firing rate.  if the average firing rate across the simulation for any neuron exceeds this value, some of that neuron's spikes will not be saved


p = paramType(PPrecision,PScale,FloatPrecision,IntPrecision,PType,seed,rng_func,example_neurons,wid,train_duration,penlambda,penlamFF,penmu,frac,learn_every,stim_on,stim_off,train_time,dt,Nsteps,Ncells,Ne,Ni,taue,taui,K,L,Lffwd,Lexc,Linh,wpscale,
je,ji,jx,mu,vre,threshe,threshi,refrac,tauedecay,tauidecay,taudecay_plastic,noise_model,sig,correlation_var,maxrate,
genStim_file, genStim_args, genTarget_file, genTarget_args, genFfwdRate_file, genFfwdRate_args, genStaticWeights_file, genStaticWeights_args, genPlasticWeights_file, genPlasticWeights_args,
choose_task_func)
