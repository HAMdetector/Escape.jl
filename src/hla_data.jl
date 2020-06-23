export HLAData

mutable struct HLAData <: AbstractHLAData
    name::String
    records::Vector{FASTX.FASTA.Record}
    hla_types::Vector{HLAType}
    phylogenetic_tree::Union{Missing, PhylogeneticTree}

    function HLAData(
        name::String, 
        records::Vector{FASTX.FASTA.Record}, 
        hla_types::Vector{HLAType},
        phylogenetic_tree::Union{Missing, PhylogeneticTree}
    )        

        if length(hla_types) != length(records)
            msg = "Vector of HLATypes has size $(length(hla_types)), " * 
                "expected $(length(records))"
            error(msg)
        end

        if !ismissing(phylogenetic_tree)
            m = Escape.matching(phylogenetic_tree, records)
            m isa Exception && throw(m)
        end
        
        new(name, records, hla_types, phylogenetic_tree)
    end
end

function HLAData(name::String, fasta_file::String, hla_types::Vector{HLAType})
    records = Escape.records(fasta_file)

    return HLAData(name, records, hla_types, missing)
end

function HLAData(
    name::String, 
    fasta_file::String, 
    hla_types::Vector{HLAType},
    tree::Union{Missing, PhylogeneticTree}
)   
    records = Escape.records(fasta_file)

    return HLAData(name, records, hla_types, tree)
end

name(data::AbstractHLAData) = getfield(data, :name)
records(data::AbstractHLAData) = getfield(data, :records)
hla_types(data::AbstractHLAData) = getfield(data, :hla_types)

function sequence_length(data::AbstractHLAData)
    records = Escape.records(data)
    sequence_length = maximum([length(FASTA.sequence(r)) for r in records])

    return sequence_length
end

function records(fasta_file::String)
    records = FASTA.Record[]

    open(FASTA.Reader, fasta_file) do reader
        for record in reader
            push!(records, record)
        end
    end

    return records
end