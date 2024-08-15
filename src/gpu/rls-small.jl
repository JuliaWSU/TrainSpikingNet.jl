function rls(itask,
             raug::CuVector{T}, k, k2, delta, Ncells, r, rX,
             Pinv, pivot, pivot64, workspace_gpu, workspace_cpu, devinfo, u_bal, utarg,
             rrXhistory, charge0, LX, penmu, penlamFF, penlambda,
             learn_seq, wpIndexIn, wpIndexConvert, wpWeightX, wpWeightIn,
             wpWeightOut, plusone, exactlyzero, PScale) where T

    for ci = 1:Ncells
        copyto_raug(raug, LX, rX, r, wpIndexIn, ci)

        generate_Pinv!(Pinv, ci, wpWeightIn, charge0, LX, penmu, penlamFF, penlambda)

        # Pinv += rrXhistory' * rrXhistory
        update_Pinv_with_rrXhistory(Pinv, view(rrXhistory.buffer,:,1:rrXhistory.nframes), LX, wpIndexIn, ci)

        # k = Pinv \ (raug / PScale)
        k .= raug ./ PScale
    	@static if p.PType == Array
            #ldiv!(k, lu!(Pinv), raug ./ PScale)
            CUSOLVER.getrf!(Pinv, pivot, devinfo, workspace_gpu)
            CUSOLVER.getrs!('N', Pinv, pivot, k, devinfo)
    	elseif p.PType == Symmetric
            #ldiv!(k, bunchkaufman!(Symmetric(Pinv)), raug / PScale)
            CUSOLVER.sytrf!('U', Pinv, pivot, devinfo, workspace_gpu)
            pivot64 .= pivot
            CUSOLVER.sytrs!('U', Pinv, pivot64, k2, devinfo, workspace_gpu, workspace_cpu)
    	end

        CUDA.@allowscalar e = view(wpWeightIn,:,ci)' * (@view raug[LX+1:end]) + u_bal[ci] - utarg[learn_seq,ci,itask]
        @static p.LX>0 && (e += view(wpWeightX,ci,:)' * rX)
        delta .= e .* k
        @static if !p.benchmark
            view(wpWeightIn,:,ci) .-= @view delta[LX+1:end]
            @static p.LX>0 && (wpWeightX[ci,:] .-= (@view delta[1:LX])')
        end
    end

    wpWeightIn2Out!(wpWeightOut, wpIndexIn, wpIndexConvert, wpWeightIn)

    push!(rrXhistory, 0)
    @static p.LX>0 && (rrXhistory[1:LX, end] = rX)
    rrXhistory[LX+1:end, end] = r

    return wpWeightIn, wpWeightOut
end

function copyto_raug(raug, LX, rX, r, wpIndexIn, ci)

    function kernel(raug, LX, rX, r, wpIndexIn, ci)
        i0 = threadIdx().x + (blockIdx().x - 1) * blockDim().x
        istride = blockDim().x * gridDim().x

        @inbounds for i=i0:istride:size(raug,1)
            raug[i] = i<=LX ? rX[i] : r[wpIndexIn[i,ci]]
        end
        return nothing
    end

    kernel = @cuda launch=false kernel(raug, LX, rX, r, wpIndexIn, ci)
    threads, blocks = configurator(kernel, size(raug,1))
    kernel(raug, LX, rX, r, wpIndexIn, ci; threads=threads, blocks=blocks)
end

function update_Pinv_with_rrXhistory(Pinv, rrXhistory, LX, wpIndexIn, ci)

    function kernel(Pinv, rrXhistory, LX, wpIndexIn, ci)
        i0 = threadIdx().x + (blockIdx().x - 1) * blockDim().x
        j0 = threadIdx().y + (blockIdx().y - 1) * blockDim().y
        istride = blockDim().x * gridDim().x
        jstride = blockDim().y * gridDim().y

        @inbounds for i=i0:istride:size(Pinv,1), j=j0:jstride:size(Pinv,2)
            irrX = i<=LX ? i : wpIndexIn[i,ci]
            jrrX = j<=LX ? j : wpIndexIn[j,ci]
            for h=1:size(rrXhistory,2)
                Pinv[i,j] += rrXhistory[irrX,h] * rrXhistory[jrrX,h]
            end
        end
        return nothing
    end

    kernel = @cuda launch=false kernel(Pinv, rrXhistory, LX, wpIndexIn, ci)
    threads, blocks = configurator(kernel, size(Pinv,1), size(Pinv,2))
    kernel(Pinv, rrXhistory, LX, wpIndexIn, ci; threads=threads, blocks=blocks)
end
