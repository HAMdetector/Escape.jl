export RateMatrix

struct RateMatrix{T <: Real} <: AbstractMatrix{T}
    m::Matrix{T}
    s::Vector{String}
    
    function RateMatrix(m::Matrix{T}, s::Vector{String}) where T <: Real
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

        new{T}(m, s)
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
Base.exp(A::RateMatrix) = RateMatrix(exp(A.m), A.s)

function transp(from_state::String, to_state::String, r::RateMatrix, branch_length::Real)
    return exp(r * branch_length)[from_state, to_state]
end

# Likelihood of observing each state at vertex k
function L_k(tree::PhylogeneticTree, k::Int, r::RateMatrix)
    states = r.s
    if isleaf(k, tree)
        state = get_property(tree, k, :state)
        state ∈ r.s || error("state $state at vertex $k found, must be one of $(r.s)")
        
        L = Dict(s => log(0) for s in states)
        L[state] = ismissing(state) ? log(0) : log(1)
        
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