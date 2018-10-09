export RateMatrix, TransitionMatrix, rate_matrix, transp, stationary

abstract type StateMatrix end

struct RateMatrix{T} <: StateMatrix
    m::Matrix{T}
    s::Vector{String}
    
    function RateMatrix(m::Matrix{T}, s::Vector{String}) where T
        length(s) == size(m)[1] || 
            throw(DimensionMismatch("matrix dimensions must match length of state names."))
        
        size(m)[1] == size(m)[2] || 
            throw(DimensionMismatch("matrix got dimensions $(size(m)), needs (2, 2)."))

        all(abs(sum(m[i,:])) < 10^-10 for i in 1:size(m)[2]) || 
            throw(ErrorException("rows of rate matrix must sum to 0."))
        
        all(diag(m) .<= 0) ||
            throw(ErrorException("diagonal matrix elements must be < 0."))

        new{T}(m, s)
    end
end

struct TransitionMatrix{T} <: StateMatrix
    m::Matrix{T}
    s::Vector{String}

    function TransitionMatrix(m::Matrix{T}, s::Vector{String}) where T
        length(s) == size(m)[1] || 
            throw(DimensionMismatch("matrix dimensions must match length of state names."))
        
        size(m)[1] == size(m)[2] || 
            throw(DimensionMismatch("matrix dimensions $(size(m)), needs (2, 2)."))

        all(abs(sum(m[i,:]) - 1) < 10^-10 for i in 1:size(m)[2]) || 
            throw(ErrorException("rows of transition matrix must sum to 1."))
        
        all(zero(T) <= e <= one(T) for e in m) ||
            throw(ErrorException("matrix elements must be between 0 and 1."))

        new{T}(m, s)
    end
end

rate_matrix(m::Matrix, s::Vector{String}) = RateMatrix(m, s)
transition_matrix(m::Matrix, s::Vector{String}) = TransitionMatrix(m, s)

function Base.getindex(A::StateMatrix, d1::String, d2::String)
    d1 in A.s && d2 in A.s || error("invalid index.")
    
    i = findfirst(x -> x == d1, A.s)
    j = findfirst(x -> x == d2, A.s)

    return A[i, j]
end

Base.size(A::StateMatrix) = size(A.m)
Base.getindex(A::StateMatrix, i::Int, j::Int) = getindex(A.m, i, j)
Base.:(*)(A::RateMatrix, x::Real) = RateMatrix(A.m * x, A.s)
Base.exp(A::RateMatrix{T}) where T = TransitionMatrix(exp(A.m), A.s)

function transp(r::RateMatrix, from_state::String, to_state::String, branch_length::Real)
    p = exp(r * branch_length)[from_state, to_state]

    return p
end

function stationary(r::RateMatrix{T}) where T
    Dict(k => v for (k,v) in zip(r.s, diag(exp(r.m * 1000))))
end