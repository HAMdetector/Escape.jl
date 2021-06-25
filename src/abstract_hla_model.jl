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