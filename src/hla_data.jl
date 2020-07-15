export HLAData

mutable struct HLAData <: AbstractHLAData
    name::String
    records::Vector{FASTX.FASTA.Record}
    hla_types::Vector{HLAType}
    tree::Union{Missing, PhylogeneticTree}
    stan_input::Union{Missing, Dict{String, Any}}

    function HLAData(
        name::String, 
        records::Vector{FASTX.FASTA.Record}, 
        hla_types::Vector{HLAType},
        tree::Union{Missing, PhylogeneticTree},
        stan_input::Union{Missing, Dict{String, Any}}
    )        

        if length(hla_types) != length(records)
            msg = "Vector of HLATypes has size $(length(hla_types)), " * 
                "expected $(length(records))"
            error(msg)
        end

        if !ismissing(tree)
            m = Escape.matching(tree, records)
            m isa Exception && throw(m)
        end
        
        if !ismissing(stan_input)
            m = ErrorException("stan_input is not in a valid format.")
            is_valid(stan_input) || throw(m)
            length(records) == stan_input["R"]
        end
        
        new(name, records, hla_types, tree, stan_input)
    end
end

function HLAData(name::String, fasta_file::String, hla_types::Vector{HLAType})
    records = Escape.records(fasta_file)

    return HLAData(name, records, hla_types, missing, missing)
end

function HLAData(
    name::String, 
    fasta_file::String, 
    hla_types::Vector{HLAType},
    tree::Union{Missing, PhylogeneticTree},
    stan_input::Union{Missing, Dict{String, Any}}
)   
    records = Escape.records(fasta_file)

    return HLAData(name, records, hla_types, tree, stan_input)
end

name(data::AbstractHLAData) = getfield(data, :name)
records(data::AbstractHLAData) = getfield(data, :records)
hla_types(data::AbstractHLAData) = getfield(data, :hla_types)
stan_input(data::AbstractHLAData) = getfield(data, :stan_input)

function is_valid(input::Dict{String, Any})
    m = ErrorException("stan_input is not in a valid format.")

    required_keys = ("Z", "rs", "S","idx", "X", "y", "N", "D", "R", "phy")
    all(map(x -> x in keys(input), required_keys)) || return false
    size(input["Z"]) == (input["R"], input["D"]) || return false
    maximum(input["rs"]) == input["R"] || return false
    maximum(input["idx"]) == input["S"] || return false
    size(input["X"]) == (input["S"], input["D"]) || return false
    size(input["phy"]) == (input["R"], input["S"]) || return false

    return true
end

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