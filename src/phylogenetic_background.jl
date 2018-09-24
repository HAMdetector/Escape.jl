export RateMatrix

struct RateMatrix# {T} <: AbstractMatrix{T}
    m#{T}
    s::Vector{String}
    
    function RateMatrix(m::Matrix, s::Vector{String})
        if size(m, 1) != size(m, 2)
            d = DimensionMismatch(string("rate matrix must have dimensions (2, 2),",
                                         " got $(size(m))"))
            throw(d)
        end

        if length(s) != size(m, 1)
            d = DimensionMismatch(string("state vector must match the number of rows",
                                         " and columns of the rate matrix."))
            throw(d)
        end

        new(m, s)
    end
end

function Base.getindex(A::RateMatrix, d1::String, d2::String)
    d1 ∈ A.s || error("$d1 is not a state in the rate matrix.")
    d2 ∈ A.s || error("$d2 is not a state in the rate matrix.")

    i = findfirst(x -> x == d1, A.s)
    j = findfirst(x -> x == d2, A.s)

    return A[i, j]
end

function Base.setindex!(A::RateMatrix, v, d1::String, d2::String)
    d1 ∈ A.s || error("$d1 is not a state in the rate matrix.")
    d2 ∈ A.s || error("$d2 is not a state in the rate matrix.")

    i = findfirst(x -> x == d1, A.s)
    j = findfirst(x -> x == d2, A.s)

    A[i, j] = v
end

Base.size(A::RateMatrix) = size(A.m)
Base.getindex(A::RateMatrix, i::Int, j::Int) = getindex(A.m, i, j)
Base.setindex!(A::RateMatrix, v, i::Int, j::Int) = setindex!(A.m, v, i, j)
Base.:(*)(A::RateMatrix, x::Real) = RateMatrix(A.m * x, A.s)
#Base.exp(A::RateMatrix) = RateMatrix(exp(A.m), A.s)

function Base.exp(A::RateMatrix)
    res = Matrix{Float64}(I, size(A)[1], size(A)[1]) + A.m

    for i in 2:20
        res .= res + A.m^i/factorial(i)
    end

    return sum(res)
end

function transp(from_state::String, to_state::String, r::RateMatrix, branch_length::Real)
    return exp(r * branch_length)[from_state, to_state]
end

# Likelihood of observing each state at vertex k
function L_k(tree::PhylogeneticTree, k::Int, r::RateMatrix)
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
    i = L_k(tree, children[1], r) # Likelihood of observing each state at child i
    j = L_k(tree, children[2], r) # Likelihood of observing each state at child j
    l_ki = get_prop(tree.graph, Edge(k => children[1]), :branch_length)
    l_kj = get_prop(tree.graph, Edge(k => children[2]), :branch_length)

    L = Dict()
    for state in states
        L_i = StatsFuns.logsumexp((log(transp(state, s, r, l_ki)) + i[s] for s in states))
        L_j = StatsFuns.logsumexp((log(transp(state, s, r, l_kj)) + j[s] for s in states))
        L[state] = L_i + L_j
    end

    return L
end

# Likelihood of observing a state pattern on a tree for given transition rates.
function L(tree::PhylogeneticTree, r::RateMatrix)
    states = r.s
    p = L_k(tree, 1, r)
    stat = exp(r * 1000)

    return StatsFuns.logsumexp(p[s] + log(stat[s, s]) for s in r.s)
end

function estimate_parameters(tree::PhylogeneticTree)
    function create_rate_matrix(x)
        x = transpose(reshape(x, sqrt(x), :))
        return RateMatrix(x, ["0", "1"])
    end
    f(x) = exp(RateMatrix(x, ["0", "1"]))
    g = grad(f)
    #r = RateMatrix([-2 2; 3 -3], ["0", "1"])
    g([-2.0 2.0; 3.0 -3.0])
    #g([-2 2; 3 -3])
end

############################################################################################
############################################################################################

function transp(from_state::String, to_state::String, branch_length::Real, π::Real, λ::Real)
    π_b = to_state == "1" ? π : 1 - π

    if from_state == to_state
        return exp(-λ*branch_length) + π_b * (1 - exp(-λ*branch_length))
    else
        return π_b*(1 - exp(-λ*branch_length))
    end
end

# Likelihood of observing each state at vertex k
function L_k(tree::PhylogeneticTree, k::Int, π::Real, λ::Real)
    states = ["0", "1"]
    if isleaf(k, tree)
        state = get_property(tree, k, :state)
        ismissing(state) && return  Dict(s => log(1) for s in states)
        state ∈ states || error("state $state at vertex $k found, must be one of $(states)")
        
        L = Dict(s => log(0) for s in states)
        L[state] = log(1)
        
        return L
    end

    children = outneighbors(tree.graph, k)
    i = L_k(tree, children[1], π, λ) # Likelihood of observing each state at child i
    j = L_k(tree, children[2], π, λ) # Likelihood of observing each state at child j
    l_ki = get_prop(tree.graph, Edge(k => children[1]), :branch_length)
    l_kj = get_prop(tree.graph, Edge(k => children[2]), :branch_length)

    L = Dict()
    for state in states
        L_i = StatsFuns.logsumexp(log(transp(state, s, l_ki, π, λ)) + i[s] for s in states)
        L_j = StatsFuns.logsumexp(log(transp(state, s, l_kj, π, λ)) + j[s] for s in states)
        L[state] = L_i + L_j
    end

    return L
end

function L(tree::PhylogeneticTree, π::Real, λ::Real)
    p = L_k(tree, 1, π, λ)
    stat = Dict("1" => π, "0" => 1 - π)

    return StatsFuns.logsumexp(p[s] + log(stat[s]) for s in ["0", "1"])
end

function parameter_estimates(tree::PhylogeneticTree)
    invlogit(x) = exp(x)/(1 + exp(x))

    f(x) = -L(tree, invlogit(x[1]), exp(x[2]))
    opt = Optim.optimize(f, [0.4, 1.1])

    return (π = invlogit(opt.minimizer[1]), λ = exp(opt.minimizer[2]))
end

function state_probabilities(tree::PhylogeneticTree)
    df = DataFrame(sequence = String[], prob = Float64[])
    ctree = deepcopy(tree)
    
    π, λ = parameter_estimates(tree)

    for leaf in leaves(tree)
        original_state = get_property(ctree, leaf, :state)

        set_property!(ctree, leaf, :state, "0")
        L_0 = L(ctree, π, λ)

        set_property!(ctree, leaf, :state, "1")
        L_1 = L(ctree, π, λ)

        prob = exp(L_1) * π / (exp(L_1) * π + exp(L_0) * (1 - π))
        push!(df, (get_property(ctree, leaf, :name), prob))

        set_property!(ctree, leaf, :state, original_state)
    end

    return df
end