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

    observed_states = map(leaf -> get_property(ctree, leaf, :state), leaves(ctree))
    count_1 = count(observed_states .=== "1")
    count_0 = count(observed_states .=== "0")
    priors = Dict{String, Float64}(
        "0" => log(count_0 / (count_0 + count_1)), 
        "1" => log(count_1 / (count_0 + count_1))
    )

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
            p[leafname][state] = exp((loglikelihoods[state]) - 
                logsumexp((loglikelihoods[s]) for s in states))
        end

        set_property!(ctree, leaf, :state, original_state)
    end

    return p
end

function L(tree::PhylogeneticTree)
    fasta_file = write_fasta(tree)
    tree_file = write_newick(tree)

    output = read(`$(raxml_ng()) --evaluate --msa $fasta_file --tree $tree_file
        --model BIN --opt-model on --opt-branches off --threads 1 --nofiles`, 
        String)

    ll = match(r"Final LogLikelihood: (.+)\n", output).captures[1]

    rm(fasta_file)
    rm(tree_file)

    return parse(Float64, ll)
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