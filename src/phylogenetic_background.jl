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
    tree = phylogenetic_tree(data)
    ctree = deepcopy(tree)
    annotate!(ctree, data, replacement)

    p = state_probabilities(ctree)
    phylogeny_effect = [p[s]["1"] for s in string.(1:length(y))]
    phylogeny_effect = [min(max(0.01, x), 0.99) for x in phylogeny_effect]

    return phylogeny_effect
end

function state_probabilities(tree::PhylogeneticTree)
    p = Dict{String, Dict{String, Float64}}()
    ctree = deepcopy(tree)

    states = ["0", "1"]

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
    fasta_file = write_fasta(tree)
    newick_file = write_newick(tree)

    try
        ll = L(fasta_file, newick_file)
        return ll
    finally
        rm(fasta_file)
        rm(newick_file)
    end
end

function L(fasta_file::String, tree_file::String)

    output = read(`raxml-ng --evaluate --msa $fasta_file --tree $tree_file 
        --model BIN --nofiles --opt-model off --opt-branches off
        --threads 1 --force`, String)
    m = match(r"Final LogLikelihood: (.+)\n", output)
    ll = tryparse(Float64, m.captures[1])
    
    !isnothing(ll) || error("raxml-ng did not return a tree loglikelihood.")

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
        for (i, s) in enumerate(states(tree))
            write(io, ">$(i)\n")
            write(io, "$s\n")
        end
    finally
        close(io)
    end

    return path
end