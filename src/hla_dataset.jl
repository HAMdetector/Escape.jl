export AbstractHLADataset, HLADataset

abstract type AbstractHLADataset end

struct HLADataset <: AbstractHLADataset
    name::String
    data::Vector{<: AbstractHLAData}
end
