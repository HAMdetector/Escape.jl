export HLAData

mutable struct FastaFile
    path::String
end

struct HLAData <: AbstractHLAData
    name::String
    fasta_file::FastaFile
    hla_types::Vector{HLAType}
    tree::Union{PhylogeneticTree, Missing}

    function HLAData(name::String, fasta_file::FastaFile, hla_types::Vector{HLAType},
            tree::Union{PhylogeneticTree, Missing})
        
        reader = BioSequences.FASTA.Reader(open(fasta_file.path, "r"))
        fasta_length = length(collect(reader))
        close(reader)

        if length(hla_types) != fasta_length
            error("Vector of HLATypes has size $(length(hla_types)), " * 
                  "expected $(fasta_length)")
        end

        if !ismissing(tree)
            matching(tree, fasta_file.path) isa Exception && 
                throw(matching(tree, fasta_file.path))
        end
        
        new(name, fasta_file, hla_types, tree)
    end
end

function HLAData(name::String, fasta_path::String, hla_types::Vector{HLAType},
                 tree::Union{PhylogeneticTree, Missing})
    fasta_file = FastaFile(fasta_path)
    
    return HLAData(name, fasta_file, hla_types, tree)
end