export AbstractHLADataset, HLADataset

abstract type AbstractHLADataset end

struct HLADataset <: AbstractHLADataset
    name::String
    data::Vector{<: AbstractHLAData}
end

function Escape.HLADataset(::Val{:Test})
    sequence_dir = joinpath(@__DIR__, "..", "data", "test_sequences")
    fasta_path = joinpath(sequence_dir, "test_large.fasta")
    
    hla_types = rand(HLAType, 15)
    hla_types[5] = HLAType((parse_allele("A11"), parse_allele("A13"),
                            parse_allele("B53"), parse_allele("B56"),
                            missing, missing))
    hla_types[6] = HLAType((missing, missing, missing, missing, missing, missing))
    hla_data_1 = HLAData("protein_1", fasta_path, hla_types, missing, missing)
    hla_data_2 = HLAData("protein_2", fasta_path, rand(HLAType, 15), missing, missing)

    if !isfile(joinpath(sequence_dir, "test_input_1.jls"))
        input = Escape.stan_input(Escape.HLAModel{4}(), hla_data_1)
        serialize(joinpath(sequence_dir, "test_input_1.jls"), input)
    end
    
    if !isfile(joinpath(sequence_dir, "test_input_2.jls"))
        input = Escape.stan_input(Escape.HLAModel{4}(), hla_data_2)
        serialize(joinpath(sequence_dir, "test_input_2.jls"), input)
    end

    input_1 = deserialize(joinpath(sequence_dir, "test_input_1.jls"))
    input_2 = deserialize(joinpath(sequence_dir, "test_input_2.jls"))

    hla_data_1 = HLAData("protein_1", fasta_path, hla_types, missing, input_1)
    hla_data_2 = HLAData("protein_2", fasta_path, rand(HLAType, 15), missing, input_2)

    dataset = HLADataset("Test", [hla_data_1, hla_data_2])

    return dataset
end