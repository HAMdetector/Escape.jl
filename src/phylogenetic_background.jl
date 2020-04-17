function phylogeny_probabilities(
    replacement::Replacement, data::AbstractHLAData, tree::PhylogeneticTree
    )

    y = targets(replacement, data)
    ctree = deepcopy(tree)
    annotate!(ctree, data, replacement)

    p = state_probabilities(ctree)
    phylogeny_effect = [p[s]["1"] for s in string.(1:length(y))]
    phylogeny_effect = [min(max(0.01, x), 0.99) for x in phylogeny_effect]

    return phylogeny_effect
end

function phylogeny_probabilities(replacement::Replacement, data::AbstractHLAData)
    y = targets(replacement, data)
    ctree = deepcopy(data.tree)
    annotate!(ctree, data, replacement)

    p = state_probabilities(ctree)
    phylogeny_effect = [p[s]["1"] for s in string.(1:length(y))]
    phylogeny_effect = [min(max(0.01, x), 0.99) for x in phylogeny_effect]

    return phylogeny_effect
end

function state_probabilities(tree::PhylogeneticTree)
    p = Dict{String, Dict{String, Float64}}()
    ctree = deepcopy(tree)

    states = [get_property(tree, x, :state) for x in leaves(ctree)] |> 
        skipmissing |> unique |> sort

    for leaf in leaves(ctree)
        leafname = get_property(ctree, leaf, :name)
        p[leafname] = Dict{String, Float64}()

        original_state = get_property(ctree, leaf, :state)
        loglikelihoods = Dict{String, Float64}()
        
        for state in states
            set_property!(ctree, leaf, :state, state)
            loglikelihoods[state] = L(ctree)
        end

        for state in states
            p[leafname][state] = exp(loglikelihoods[state] - 
                logsumexp(loglikelihoods[s] for s in states))
        end

        set_property!(ctree, leaf, :state, original_state)
    end

    return p
end

function L(tree::PhylogeneticTree)
    fasta_filepath = write_fasta(tree)
    newick_filepath = write_newick(tree)

    try
        ll = L(fasta_filepath, newick_filepath)
        return ll
    finally
        rm(fasta_filepath)
        rm(newick_filepath)
    end
end

function L(fasta_filepath::String, tree_filepath::String)

    output = read(`raxml-ng --evaluate --msa $fasta_filepath --tree $tree_filepath 
        --model BIN --nofiles --opt-model off --opt-branches off
        --threads 1`, String)
    m = match(r"Final LogLikelihood: (.+)\n", output)
    ll = tryparse(Float64, m.captures[1])
    
    !isnothing(ll) || throw(error("raxml-ng did not return a tree loglikelihood."))

    return ll
end

function write_newick(tree::PhylogeneticTree)
    path, io = mktemp(cleanup = true)
    try
        write(io, newick_string(tree))
    finally
        close(io)
    end

    return path
end

function write_fasta(tree::PhylogeneticTree)
    path, io = mktemp(cleanup = true)

    try
        writer = FASTA.Writer(io)

        for (i, s) in enumerate(states(tree))
            write(io, ">$(i)\n")
            write(io, "$s\n")
        end
    finally
        close(io)
    end

    return path
end