export AbstractHLAData, HLAData, AbstractHLADataset, HLADataset

abstract type AbstractHLAData end
abstract type AbstractHLADataset end

struct HLAData <: AbstractHLAData
    name::String
    fasta_file::String
    hla_types::Vector{HLAType}

    function HLAData(name::String, fasta_file::String, hla_types::Vector{HLAType})
        reader = BioSequences.FASTA.Reader(open(fasta_file, "r"))
        fasta_length = length(collect(reader))
        close(reader)

        if length(hla_types) != fasta_length
            error("Vector of HLATypes has size $(length(hla_types)), " * 
                  "expected $(fasta_length)")
        end

        new(name, fasta_file, hla_types)
    end
end 

struct HLADataset <: AbstractHLADataset
    name::String
    data::Vector{<: AbstractHLAData}
end

function hla_matrix(data::Vector{HLAType})
    alleles = limit_hla_accuracy.(sort(unique_alleles(data)))
    m = zeros(Int, length(data), length(alleles))

    for i in eachindex(data)
        hla_type = data[i]
        for allele in hla_type.alleles
            ismissing(allele) && continue
            allele = limit_hla_accuracy(allele)
            m[i, findfirst(x -> x == allele, alleles)] += 1
        end
    end

    return m
end
