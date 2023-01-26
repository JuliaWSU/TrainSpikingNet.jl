using Revise

using SparseArrays
using ProgressMeter

function potjans_params()
    conn_probs = [[0.1009,  0.1689, 0.0437, 0.0818, 0.0323, 0.,     0.0076, 0.    ],
                [0.1346,   0.1371, 0.0316, 0.0515, 0.0755, 0.,     0.0042, 0.    ],
                [0.0077,   0.0059, 0.0497, 0.135,  0.0067, 0.0003, 0.0453, 0.    ],
                [0.0691,   0.0029, 0.0794, 0.1597, 0.0033, 0.,     0.1057, 0.    ],
                [0.1004,   0.0622, 0.0505, 0.0057, 0.0831, 0.3726, 0.0204, 0.    ],
                [0.0548,   0.0269, 0.0257, 0.0022, 0.06,   0.3158, 0.0086, 0.    ],
                [0.0156,   0.0066, 0.0211, 0.0166, 0.0572, 0.0197, 0.0396, 0.2252],
                [0.0364,   0.001,  0.0034, 0.0005, 0.0277, 0.008,  0.0658, 0.1443]]

    columns_conn_probs = [col for col in eachcol(conn_probs)][1]
    layer_names = ["23E","23I","4E","4I","5E", "5I", "6E", "6I"]

    ccuf = Dict(
        k=>v for (k,v) in zip(layer_names,columns_conn_probs)
    )

    ccu = Dict("23E"=>20683, "23I"=>5834,
                "4E"=>21915, "4I"=>5479,
                "5E"=>4850, "5I"=>1065,
                "6E"=>14395, "6I"=>2948)

    ccu = Dict((k,ceil(Int64,v/35.0)) for (k,v) in pairs(ccu))
    cumulative = Dict() 
    v_old=1
    for (k,v) in pairs(ccu)
        cumulative[k]=collect(v_old:v+v_old)
        v_old=v+v_old
    end
    return (cumulative,ccu,ccuf,layer_names,columns_conn_probs,conn_probs)
end


function re_write_weights(Ncells)
    (cumulative,ccu,ccuf,layer_names,columns_conn_probs,conn_probs) = potjans_params()    
    w_mean = 87.8e-3  # nA
    ###
    # Efficiency right!
    # Lower memory footprint motivations.
    # a sparse matrix can be stored as a smaller dense matrix.
    # A 2D matrix should be stored as 1D matrix of srcs,tgts
    # A weight matrix should be stored as 1 matrix, which is redistributed in loops using 
    # the 1D matrix of srcs,tgts.
    ###
    Ncells = sum([i for i in values(ccu)])+1#max([max(i[:]) for i in values(cumulative)])
    

    w0Index_ = spzeros(Int,Ncells,Ncells)
    w0Weights = spzeros(Float32,Ncells,Ncells)
    edge_dict = Dict() # Int64,Array{Int64}
    polarity = Dict()

    for src in 1:p.Ncells
        edge_dict[src] = Int64[]
        polarity[src] = ""
    end
    Ne = 0 
    Ni = 0
    @showprogress for (i,(k,v)) in enumerate(pairs(cumulative))
        for src in v
            for (j,(k1,v1)) in enumerate(pairs(cumulative))
                for tgt in v1
                    if src!=tgt
                        prob = conn_probs[i][j]
                        if rand()<prob
                            if occursin("E",k) 
                                if occursin("E",k1)                    
                                    w0Weights[tgt,src] = p.jee#*#350.0#)/2.0
                                elseif occursin("I",k1)                    
                                    w0Weights[tgt,src] = p.jei#*150.0  
                                end
                                polarity[src]="E"
                    
                            elseif occursin("I",k)
                                if occursin("E",k1)                    
                                    w0Weights[tgt,src] = -p.jie  
                                elseif occursin("I",k1)                    
                                    w0Weights[tgt,src] = -p.jii#*2.0  
                                end
                                    polarity[src]="I"

                            end
                            append!(edge_dict[src],tgt)
                            w0Index_[tgt,src] = tgt

                        end
                    end
                end
            end
        end
    
    end
    for val in values(polarity)
        if val=="E"
            Ne+=1
        elseif val=="I"
            Ni+=1
        end
    end
    return (edge_dict,w0Weights,w0Index_,Ne,Ni)
end

function convert_dense_matrices(p)
    (edge_dict,w0Weights,w0Index_,Ne,Ni) = re_write_weights(p.Ncells)
    
    ##
    # nc0Max is the maximum number of post synaptic targets
    # its a limit on the outdegree.
    # if this is not known upfront it can be calculated on the a pre-exisiting adjacency matrix as I do below.
    ##

    nc0Max = 1

    for (k,v) in pairs(edge_dict)
        templength = length(v)
        if templength>nc0Max
            nc0Max=templength
        end
    end

    #nc0Max = Ncells-1 # outdegree
    nc0 = Int.(nc0Max*ones(Ncells))
    w0Index = spzeros(Int,nc0Max,Ncells)
    for pre_cell = 1:p.Ncells
        post_cells = edge_dict[pre_cell]
        w0Index[1:length(edge_dict[pre_cell]),pre_cell] = post_cells
    end

    return w0Index, w0Weights, nc0, w0Index_, Ne, Ni
end

function genWeights(p)

    nc0Max = Int(p.Ncells*p.pree) # outdegree
    nc0 = Int.(nc0Max*ones(p.Ncells))
    w0Index = spzeros(Int,nc0Max,p.Ncells)
    w0Weights = spzeros(nc0Max,p.Ncells)
    for precell = 1:p.Ncells
        postcells = filter(x->x!=i, collect(1:p.Ncells)) # remove autapse
        ###
        # take a small random subset of all possible post synaptic cells
        # by randomly shuffling all of the cell indexs, and then taking a small finite slice of the subset.
        ###
        tgts = sort(shuffle(postcells)[1:nc0Max]) # fixed outdegree nc0Max
        w0Index[1:nc0Max,precell] = tgts
        nexc = sum(w0Index[1:nc0Max,i] .<= p.Ne) # number of exc synapses
        if i <= p.Ne
            w0Weights[1:nexc,precell] .= p.jee  ## EE weights
            w0Weights[nexc+1:nc0Max,i] .= p.jie  ## IE weights
        else
            w0Weights[1:nexc,precell] .= p.jei  ## EI weights
            w0Weights[nexc+1:nc0Max,i] .= p.jii  ## II weights
        end
    end

    return w0Index, w0Weights, nc0
end

    
    

function genWeights_square(p)
    w_mean = 87.8e-3  # nA

    nc0Max = Int(p.Ncells*p.pree) # outdegree
    nc0 = Int.(p.Ncells*ones(p.Ncells))
    w0Index = spzeros(Int,p.Ncells,p.Ncells)
    w0Weights = spzeros(p.Ncells,p.Ncells)
    for i = 1:p.Ncells
        postcells = filter(x->x!=i, collect(1:p.Ncells)) # remove autapse
        w0Index[1:nc0Max,i] = sort(shuffle(postcells)[1:nc0Max]) # fixed outdegree nc0Max
        nexc = sum(w0Index[1:nc0Max,i] .<= p.Ne) # number of exc synapses
        if i <= p.Ne
            w0Weights[1:nexc,i] .= w_mean  ## EE weights
            w0Weights[nexc+1:nc0Max,i] .= p.jie  ## IE weights
        else
            w0Weights[1:nexc,i] .= p.jei  ## EI weights
            w0Weights[nexc+1:nc0Max,i] .= p.jii  ## II weights
        end
    end
    
    return w0Index, w0Weights, nc0
    
    
    end

