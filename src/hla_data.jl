export HLAData

mutable struct HLAData <: AbstractHLAData
    name::String
    fasta_file::String
    hla_types::Vector{HLAType}
    tree::Union{Missing, PhylogeneticTree}

    function HLAData(
        name::String, 
        fasta_file::String, 
        hla_types::Vector{HLAType},
        tree::Union{PhylogeneticTree, Missing}
    )        
        open(FASTA.Reader, fasta_file) do reader
            fasta_length = length(collect(reader))

            if length(hla_types) != fasta_length
                error("Vector of HLATypes has size $(length(hla_types)), " * 
                      "expected $(fasta_length)")
            end
        end

        if !ismissing(tree) && matching(tree, fasta_file) isa Exception
            throw(matching(tree, fasta_file))
        end
        
        new(name, fasta_file, hla_types, tree)
    end
end

function HLAData(name::String, fasta_path::String, hla_types::Vector{HLAType})
    return HLAData(name, fasta_path, hla_types, missing)
end