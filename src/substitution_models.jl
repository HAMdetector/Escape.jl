export SubstitutionModel, TwoStateGTR

abstract type SubstitutionModel end

struct TwoStateGTR <: SubstitutionModel
    α::Float64
    π_1::Float64
    π_2::Float64

    function TwoStateGTR(α, π_1, π_2)
        1 - abs(π_1 + π_2) < 0.0001 || 
            throw(ErrorException("π_1 + π_2 must be 1, got $(π_1 + π_2).")) 
        new(α, π_1, π_2)
    end
end

TwoStateGTR(; α = nothing, π_1 = nothing, π_2 = nothing) = TwoStateGTR(α, π_1, π_2)

function rate_matrix(model::TwoStateGTR, states::Vector{String})
    α = model.α
    π_1 = model.π_1
    π_2 = model.π_2

    m = [-α*π_1 α*π_1 ; α*π_2 -α*π_2]
    rate_matrix(m, states)
end

