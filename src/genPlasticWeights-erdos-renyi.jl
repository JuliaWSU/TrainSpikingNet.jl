#=
return three matrices and a vector specifying the plastic connectivity.

the first returned matrix (called wpWeightFfwd in the source code) is
Ncells x Lffwd and contains the (initial) weights of the feed forward
presynaptic neurons.

the second (wpWeightIn, Ncells columns) and third (wpIndexIn, Ncells rows)
matrices contain the weights and indices of the recurrent presynaptic
neurons.  the final returned variable is a vector of length Ncells (ncpIn)
which specifies how many presynaptic connections each neuron has.
=#

function genPlasticWeights(args, w0Index, nc0, ns0)
    Ncells, frac, Ne, L, Lexc, Linh, Lffwd, wpee, wpie, wpei, wpii, wpffwd = map(x->args[x],
            [:Ncells, :frac, :Ne, :L, :Lexc, :Linh, :Lffwd, :wpee, :wpie, :wpei, :wpii, :wpffwd])

    # order neurons by their firing rate
    frac_neurons_selected = p.frac
    frac_cells = round(Int, frac_neurons_selected*p.Ne)
    exc_ns0 = ns0[1:p.Ne]
    inh_ns0 = ns0[p.Ne+1:p.Ncells]
    exc_ordered = sortperm(exc_ns0)
    inh_ordered = collect(p.Ne+1:p.Ncells)[sortperm(inh_ns0)]
    exc_selected = sort(exc_ordered[end-frac_cells+1:end])
    inh_selected = sort(inh_ordered[end-frac_cells+1:end])
    
    # define weights_plastic
    wpWeightIn = Array{Float64}(undef, p.Lexc+p.Linh, p.Ncells)
    wpIndexIn = Array{Int}(undef, p.Ncells, p.Lexc+p.Linh)
    ncpIn = Array{Int}(undef, p.Ncells)
    for postCell = 1:p.Ncells
        # select random exc and inh presynaptic neurons
        # (1) select consecutive neurons
        # rnd_start = rand(1:length(exc_selected)-p.L+1)
        # indE = sort(exc_selected[rnd_start:rnd_start+p.L-1])
        # indI = sort(inh_selected[rnd_start:rnd_start+p.L-1])

        # (2) select random neurons
        indE = sample(rng, exc_selected, p.L, replace=false, ordered=true)
        indI = sample(rng, inh_selected, p.L, replace=false, ordered=true)

        # build wpIndexIn
        ind  = [indE; indI]
        wpIndexIn[postCell,:] = ind
        ncpIn[postCell] = length(ind)

        # initial exc and inh plastic weights
        if postCell <= p.Ne
            wpWeightIn[1:p.Lexc, postCell] .= wpee
            wpWeightIn[p.Lexc.+(1:p.Linh), postCell] .= wpei
        else
            wpWeightIn[1:p.Lexc, postCell] .= wpie
            wpWeightIn[p.Lexc.+(1:p.Linh), postCell] .= wpii
        end
    end

    # define feedforward weights to all neurons
    #       - wpWeightFfwd = randn(p.Ncells, p.Lffwd) * wpffwd
    #       - initial weights, wpffwd = 0
    wpWeightFfwd = randn(rng, p.Ncells, p.Lffwd) * wpffwd
    
    return wpWeightFfwd, wpWeightIn, wpIndexIn, ncpIn
end
