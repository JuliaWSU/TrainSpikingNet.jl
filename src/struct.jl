struct paramType
    PPrecision::DataType
    PScale::Int
    FloatPrecision::DataType
    IntPrecision::DataType
    PType::UnionAll  #
    seed::Union{Nothing,Int}
    rng_func::Dict
    example_neurons::Int
    wid::Int
    train_duration::Float64
    penlambda::Float64
    penlamFF::Float64
    penmu::Float64
    frac::Float64
    learn_every::Float64
    stim_on::Float64
    stim_off::Float64
    train_time::Float64
    dt::Float64
    Nsteps::Int64
    Ncells::Int64
    Ne::Int64
    Ni::Int64
    taue::Float64
    taui::Float64
    K::Int64
    L::Int64
    Lffwd::Int64
    Lexc::Int64
    Linh::Int64
    wpscale::Float64
    je::Float64
    ji::Float64
    jx::Float64
    mu::Vector{Float64}
    vre::Float64
    threshe::Float64
    threshi::Float64
    refrac::Float64
    tauedecay::Float64
    tauidecay::Float64
    taudecay_plastic  #::Union{Float64,Vector{Float64}}
    noise_model::Symbol
    sig::Float64
    correlation_var::Symbol
    maxrate::Float64
    genStim_file::AbstractString
    genStim_args::AbstractDict
    genTarget_file::AbstractString
    genTarget_args::AbstractDict
    genFfwdRate_file::AbstractString
    genFfwdRate_args::AbstractDict
    genStaticWeights_file::AbstractString
    genStaticWeights_args::AbstractDict
    genPlasticWeights_file::AbstractString
    genPlasticWeights_args::AbstractDict
    choose_task_func::Expr
end
