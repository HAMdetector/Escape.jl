export PhylogenyHLAData

struct PhylogenyHLAData <: AbstractHLAData
    name::String
    fasta_file::String
    hla_types::Vector{HLAType}
    tree::PhylogeneticTree

    function PhylogenyHLAData(name::String, fasta_file::String, hla_types::Vector{HLAType},
                              tree::PhylogeneticTree)
        reader = BioSequences.FASTA.Reader(open(fasta_file, "r"))
        fasta_length = length(collect(reader))
        close(reader)

        if length(hla_types) != fasta_length
            error("Vector of HLATypes has size $(length(hla_types)), " * 
                    "expected $(fasta_length)")
        end

        if length(leaves(tree)) != fasta_length
            error("Phylogenetic tree has $(lenght(leaves(tree))) leaves, " * 
                "expected $(fasta_length)")
        end

        new(name, fasta_file, hla_types, tree) 
    end
end

function phylogenetic_tree(data::PhylogenyHLAData)
    return data.tree
end