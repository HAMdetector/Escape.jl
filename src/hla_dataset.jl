export AbstractHLADataset, HLADataset

abstract type AbstractHLADataset end

struct HLADataset <: AbstractHLADataset
    name::String
    data::Vector{<: AbstractHLAData}
end

function Escape.HLADataset(::Val{:Test})
    fasta_path = joinpath(@__DIR__, "..", "data", "test_sequences", "test_large.fasta")
    hla_types = rand(HLAType, 15)
    hla_types[5] = HLAType((parse_allele("A11"), parse_allele("A13"),
                            parse_allele("B53"), parse_allele("B56"),
                            missing, missing))
    hla_types[6] = HLAType((missing, missing, missing, missing, missing, missing))
    hla_data_1 = HLAData("protein_1", fasta_path, hla_types, missing, missing)
    hla_data_2 = HLAData("protein_2", fasta_path, rand(HLAType, 15), missing, missing)

    dataset = HLADataset("Test", [hla_data_1, hla_data_2])

    return dataset
end