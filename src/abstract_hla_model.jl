export AbstractHLAModel
export AbstractHLAModelResult
export HLAModel
export HLAModelResult

abstract type AbstractHLAModel end
abstract type AbstractHLAModelResult end

struct HLAModel{x} <: AbstractHLAModel end

struct HLAModelResult{T}
    model::HLAModel{T}
    data::AbstractHLAData
    sf::StanInterface.Stanfit
    replacements::Vector{Replacement}
    alleles::Vector{HLAAllele}
end

struct PhylogenyIncluded end
struct PhylogenyExluded end

struct EpitopesIncluded end
struct EpitopesExcluded end

phylogeny_information(::HLAModel) = PhylogenyIncluded()
epitope_information(::HLAModel) = EpitopesIncluded()