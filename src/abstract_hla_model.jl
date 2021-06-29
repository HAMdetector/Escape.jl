export AbstractHLAModel
export HLAModel
export HLAModelResult

abstract type AbstractHLAModel end

struct HLAModel{x} <: AbstractHLAModel end

struct HLAModelResult{T}
    model::HLAModel{T}
    data::AbstractHLAData
    sf::StanInterface.Stanfit
    replacements::Vector{Replacement}
    alleles::Vector{HLAAllele}
end

function Base.show(io::IO, x::HLAModelResult{T}) where T
    print(io, 
        "HLAModelResult (model $T): $(name((x.data)))\n" *
        "number of replacements: $(length(replacements(x)))\n" *
        "number of alleles:      $(length(x.alleles))"
    )
end