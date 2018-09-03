export AbstractHLAData, HLAData

abstract type AbstractHLAData end

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