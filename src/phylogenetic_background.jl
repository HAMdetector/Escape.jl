function phylogeny_probabilities(
    replacement::Replacement, data::AbstractHLAData, tree::PhylogeneticTree
    )

    y = targets(replacement, data)
    ctree = deepcopy(tree)
    annotate!(ctree, data, replacement)

    p = state_probabilities(ctree, TwoStateGTR)
    phylogeny_effect = [p[s]["1"] for s in string.(1:length(y))]
    phylogeny_effect = [min(max(0.01, x), 0.99) for x in phylogeny_effect]

    return phylogeny_effect
end

function phylogeny_probabilities(replacement::Replacement, data::AbstractHLAData)
    y = targets(replacement, data)
    ctree = deepcopy(data.tree)
    annotate!(ctree, data, replacement)

    p = state_probabilities(ctree, TwoStateGTR)
    phylogeny_effect = [p[s]["1"] for s in string.(1:length(y))]
    phylogeny_effect = [min(max(0.01, x), 0.99) for x in phylogeny_effect]

    return phylogeny_effect
end

function state_probabilities(tree::PhylogeneticTree, model::Type{TwoStateGTR})
    p = Dict{String, Dict{String, Float64}}()
    r = infer_rate_matrix(tree, model)
    stat = stationary(r)
    ctree = deepcopy(tree)

    states = [get_property(tree, x, :state) for x in leaves(ctree)] |> 
        skipmissing |> unique |> sort

    for leaf in leaves(ctree)
        leafname = get_property(ctree, leaf, :name)
        p[leafname] = Dict{String, Float64}()

        original_state = get_property(ctree, leaf, :state)
        likelihoods = Dict{String, Float64}()
        
        for state in states
            set_property!(ctree, leaf, :state, state)
            likelihoods[state] = exp(L(ctree, r)) * stat[state]
        end

        for state in states
            p[leafname][state] = likelihoods[state] / sum(likelihoods[s] for s in states)
        end

        set_property!(ctree, leaf, :state, original_state)
    end

    return p
end

function infer_rate_matrix(tree::PhylogeneticTree, model::Type{TwoStateGTR})
    function l(α, π_1, π_2)
        r = rate_matrix(TwoStateGTR(α = α, π_1 = π_1, π_2 = π_2), states)
        return L(tree, r)
    end

    function ∇l(g, α, π_1, π_2)
        r = rate_matrix(TwoStateGTR(α = α, π_1 = π_1, π_2 = π_2), states)
        f = x -> l(x[1], x[2], x[3])

        g .= Calculus.gradient(f, [α, π_1, π_2])
    end

    states = ["0", "1"]

    m = Model(optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0))
    JuMP.register(m, :l, 3, l, ∇l)
    @variable(m, α, start = 1)
    @variable(m, π_1, start = 0.5)
    @variable(m, π_2, start = 0.5)

    @constraint(m, π_1 >= 0.001)
    @constraint(m, π_2 >= 0.001)
    @constraint(m, α >= 0.001)
    @constraint(m, π_1 + π_2 == 1)
    @NLconstraint(m, 2*α*π_1*π_2 == 1)

    @NLobjective(m, Max, l(α, π_1, π_2))
    optimize!(m)

    estimate = TwoStateGTR(α = JuMP.value(α), π_1 = JuMP.value(π_1), π_2 = JuMP.value(π_2))
    
    return rate_matrix(estimate, states)
end

# Likelihood of observing a state pattern on a tree for given transition rates.
function L(tree::PhylogeneticTree, r::RateMatrix)
    states = r.s
    p = L_k(tree, r, 1)
    stat = Dict(k => v for (k,v) in zip(r.s, diag(exp(r.m * 1000))))

    return StatsFuns.logsumexp(p[s] + log(stat[s]) for s in r.s)
end

# Likelihood of subtree of node k for each state
function L_k(tree::PhylogeneticTree, r::RateMatrix, k::Int)
    states = r.s
    if isleaf(k, tree)
        state = get_property(tree, k, :state)
        ismissing(state) && return  Dict(s => log(1) for s in states)
        state ∈ r.s || error("state $state at vertex $k found, must be one of $(r.s)")
        
        L = Dict(s => log(0) for s in states)
        L[state] = log(1)
        
        return L
    end

    children = outneighbors(tree.graph, k)
    i = L_k(tree, r, children[1]) # Likelihood of observing each state at child i
    j = L_k(tree, r, children[2]) # Likelihood of observing each state at child j
    l_ki = get_prop(tree.graph, Edge(k => children[1]), :branch_length)
    l_kj = get_prop(tree.graph, Edge(k => children[2]), :branch_length)

    L = Dict()
    for state in states
        L_i = StatsFuns.logsumexp(log(transp(r, state, s, l_ki)) + i[s] for s in states)
        L_j = StatsFuns.logsumexp(log(transp(r, state, s, l_kj)) + j[s] for s in states)
        L[state] = L_i + L_j
    end

    return L
end