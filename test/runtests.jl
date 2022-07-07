using Test, JLD2, SymmetricFormats

function compare_cpu_to_gpu(kind)
    init_out = readlines(`$(Base.julia_cmd()) -t 2
                          $(joinpath(@__DIR__,"..","src","init.jl"))
                          $(joinpath(@__DIR__,"cpu-$kind"))`)
    for iHz in findall(contains("Hz"), init_out)
        @test 0 < parse(Float64, match(r"([.0-9]+) Hz", init_out[iHz]).captures[1]) < 15
    end

    cp(joinpath(@__DIR__,"cpu-$kind"), joinpath(@__DIR__,"gpu-$kind"))
    cpu_out = readlines(`$(Base.julia_cmd()) -t 2
                         $(joinpath(@__DIR__,"..","src","cpu","train.jl"))
                         $(joinpath(@__DIR__,"cpu-$kind"))`)
    gpu_out = readlines(`$(Base.julia_cmd())
                         $(joinpath(@__DIR__,"..","src","gpu","train.jl"))
                         $(joinpath(@__DIR__,"gpu-$kind"))`)

    iHz = findlast(contains("Hz"), cpu_out)
    @test 0 < parse(Float64, match(r"([.0-9]+) Hz", cpu_out[iHz]).captures[1]) < 15
    iHz = findlast(contains("Hz"), gpu_out)
    @test 0 < parse(Float64, match(r"([.0-9]+) Hz", gpu_out[iHz]).captures[1]) < 15

    cpu_P = load(joinpath(@__DIR__, "cpu-$kind", "P-ckpt1.jld2"), "P")
    gpu_P = load(joinpath(@__DIR__, "gpu-$kind", "P-ckpt1.jld2"), "P")
    @test isapprox(kind=="SymmetricPacked" ? cat((x.tri for x in cpu_P)..., dims=2) :
                                             cat(cpu_P..., dims=3),
                   gpu_P)

    cpu_wpWeightIn = load(joinpath(@__DIR__, "cpu-$kind", "wpWeightIn-ckpt1.jld2"),
                          "wpWeightIn")
    gpu_wpWeightIn = load(joinpath(@__DIR__, "gpu-$kind", "wpWeightIn-ckpt1.jld2"),
                          "wpWeightIn")
    @test isapprox(cpu_wpWeightIn, gpu_wpWeightIn)
end

@testset "$kind" for kind in ("Array", "Symmetric", "SymmetricPacked")
    mkdir(joinpath(@__DIR__, "cpu-$kind"))
    open(joinpath(@__DIR__, "cpu-$kind", "param.jl"), "w") do fileout 
        for line in readlines(joinpath(@__DIR__,"param.jl"))
            if startswith(line, "PType")
                println(fileout, "PType=$kind")
            else
                println(fileout, line)
            end
        end
    end

    compare_cpu_to_gpu(kind)

    if kind!="Array"
        cpu_wpWeightIn_Array = load(joinpath(@__DIR__, "cpu-Array", "wpWeightIn-ckpt1.jld2"),
                                    "wpWeightIn")
        cpu_wpWeightIn_kind = load(joinpath(@__DIR__, "cpu-$kind", "wpWeightIn-ckpt1.jld2"),
                                   "wpWeightIn")
        @test isapprox(cpu_wpWeightIn_Array, cpu_wpWeightIn_kind)
    end
end

@testset "test" begin
    run(pipeline(`$(Base.julia_cmd())
                  $(joinpath(@__DIR__,"..","src","cpu","test.jl"))
                  $(joinpath(@__DIR__,"cpu-Array"))`, stdout=devnull))

    run(pipeline(`$(Base.julia_cmd())
                  $(joinpath(@__DIR__,"..","src","gpu","test.jl"))
                  $(joinpath(@__DIR__,"gpu-Array"))`, stdout=devnull))

    dcpu = load(joinpath(@__DIR__,"cpu-Array/test.jld2"))
    dgpu = load(joinpath(@__DIR__,"gpu-Array/test.jld2"))
    @test dcpu["nss"] == dgpu["nss"]
    @test dcpu["ineurons_to_plot"] == dgpu["ineurons_to_plot"]
    @test isapprox(dcpu["xtotals"], dgpu["xtotals"])
    @test isapprox(dcpu["timess"], dgpu["timess"])
end

@testset "pree=$pree, sig0=$sig0" for pree in [0.1, 0.0], sig0 in [0.65, 0.0]
    kind = string("pree", pree, "-sig0", sig0)
    mkdir(joinpath(@__DIR__, "cpu-$kind"))
    open(joinpath(@__DIR__, "cpu-$kind", "param.jl"), "w") do fileout 
        for line in readlines(joinpath(@__DIR__,"param.jl"))
            if contains(line, ":pree => 0.1")
                println(fileout, replace(line, ":pree => 0.1" => ":pree => $pree"))
            elseif startswith(line, "sig0 =")
                println(fileout, "sig0 = $sig0")
            elseif startswith(line, "L = ") && pree==0.0
                println(fileout, "L = 14")
            else
                println(fileout, line)
            end
        end
    end
    compare_cpu_to_gpu(kind)
end

@testset "voltage noise model" begin
    mkdir(joinpath(@__DIR__, "cpu-voltagenoise"))
    open(joinpath(@__DIR__, "cpu-voltagenoise", "param.jl"), "w") do fileout 
        for line in readlines(joinpath(@__DIR__,"param.jl"))
            if contains(line, "noise_model=:current")
                println(fileout, replace(line, "noise_model=:current" => "noise_model=:voltage"))
            else
                println(fileout, line)
            end
        end
    end
    compare_cpu_to_gpu("voltagenoise")
end

@testset "Ricciardi" begin
    mkdir(joinpath(@__DIR__, "ricciardi"))
    open(joinpath(@__DIR__, "ricciardi", "param.jl"), "w") do fileout 
        for line in readlines(joinpath(@__DIR__,"param.jl"))
            if startswith(line, "Ncells")
                println(fileout, "Ncells=4096")
            else
                println(fileout, line)
            end
        end
    end
    psth = Float64[1:100 100:-1:1 vcat(1:50,50:-1:1) (25 .+ 25*sin.(range(0,2*pi,100)))]
    save(joinpath(@__DIR__, "ricciardi", "spikerates.jld2"), "psth", psth)
    init_out = readlines(pipeline(`$(Base.julia_cmd()) -t 2
                                   $(joinpath(@__DIR__,"..","src","init.jl"))
                                   -s $(joinpath(@__DIR__,"ricciardi","spikerates.jld2"))
                                   $(joinpath(@__DIR__,"ricciardi"))`))
    xtarg = load(joinpath(@__DIR__,"ricciardi","xtarg.jld2"))["xtarg"]
    dx = diff(xtarg, dims=1)
    @test all(dx[:,1].>0)
    @test all(dx[:,2].<0)
    @test all(dx[1:49,3].>0)
    @test all(dx[51:99,3].<0)
    @test all(dx[1:25,4].>0)
    @test all(dx[26:74,4].<0)
    @test all(dx[75:99,4].>0)
end

@testset "Int16" begin
    mkdir(joinpath(@__DIR__, "Int16"))
    open(joinpath(@__DIR__, "Int16", "param.jl"), "w") do fileout 
        for line in readlines(joinpath(@__DIR__,"param.jl"))
            if startswith(line, "PPrecision")
                println(fileout, "PPrecision=Int16")
            elseif startswith(line, "PScale")
                println(fileout, "PScale=2^14")
            else
                println(fileout, line)
            end
        end
    end

    init_out = readlines(pipeline(`$(Base.julia_cmd()) -t 2
                                   $(joinpath(@__DIR__,"..","src","init.jl"))
                                   $(joinpath(@__DIR__,"Int16"))`))
    for iHz in findall(contains("Hz"), init_out)
      @test 0 < parse(Float64, match(r"([.0-9]+) Hz", init_out[iHz]).captures[1]) < 10
    end

    gpu_out = readlines(pipeline(`$(Base.julia_cmd())
                                  $(joinpath(@__DIR__,"..","src","gpu","train.jl"))
                                  $(joinpath(@__DIR__,"Int16"))`))
    iHz = findlast(contains("Hz"), gpu_out)
    @test 0 < parse(Float64, match(r"([.0-9]+) Hz", gpu_out[iHz]).captures[1]) < 10
end

@testset "feed forward" begin
    mkdir(joinpath(@__DIR__, "cpu-Lffwd"))
    open(joinpath(@__DIR__, "cpu-Lffwd", "param.jl"), "w") do fileout 
        for line in readlines(joinpath(@__DIR__,"param.jl"))
            if startswith(line, "Lffwd")
                println(fileout, "Lffwd = L>>1")
            else
                println(fileout, line)
            end
        end
    end
    compare_cpu_to_gpu("Lffwd")
end
