export HLAData

mutable struct HLAData <: AbstractHLAData
    name::String
    records::Vector{FASTX.FASTA.Record}
    hla_types::Vector{HLAType}
    tree::Union{Missing, PhylogeneticTree}
    stan_input::Union{Missing, Dict{String, Any}}
    allele_depth::Int
    replacement_mincount::Int

    function HLAData(
        name::String, 
        records::Vector{FASTX.FASTA.Record}, 
        hla_types::Vector{HLAType},
        tree::Union{Missing, PhylogeneticTree},
        stan_input::Union{Missing, Dict{String, Any}},
        allele_depth::Int,
        replacement_mincount::Int
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
        end
        
        if !(allele_depth == 1 || allele_depth == 2)
            throw(ErrorException("allele_depth must be 1 (for 2-digit accuracy " *
                "or 2 (for 4-digit accuravy)."))
        end

        if !(replacement_mincount > 0)
            throw(ErrorException("replacement_mincount must be at least 1."))
        end

        if replacement_mincount > length(records)
            throw(ErrorException("replacement_mincount must not exceed the number of" *
                "available sequences."))
        end

        new(name, records, hla_types, tree, stan_input, allele_depth, replacement_mincount)
    end
end

function HLAData(name::String, fasta_file::String, hla_types::Vector{HLAType};
        allele_depth::Int = 1, replacement_mincount::Int = 1)
    
    records = Escape.records(fasta_file)

    return HLAData(name, records, hla_types, missing, missing, allele_depth,
        replacement_mincount)
end

function HLAData(
    name::String, 
    fasta_file::String, 
    hla_types::Vector{HLAType},
    tree::Union{Missing, PhylogeneticTree},
    stan_input::Union{Missing, Dict{String, Any}};
    allele_depth::Int = 1,
    replacement_mincount::Int = 1,
)   
    records = Escape.records(fasta_file)

    return HLAData(name, records, hla_types, tree, stan_input, allele_depth, 
        replacement_mincount)
end

function HLAData(;
    alignment_file::String = "",
    tree_file::String = "",
    hla_annotation_file::String = "",
    name = splitext(basename(alignment_file))[1],
    allele_depth = 1,
    replacement_mincount = 1
)

    tree = phylogenetic_tree(readline(tree_file))
    if hla_annotation_file == ""
        check_inputs(alignment_file, tree_file)
        hla_types = let 
            r = records(alignment_file)
            hla_types = HLAType[]

            for record in r
                annotation = FASTA.identifier(record)
                push!(hla_types, HLAType(Tuple(parse_allele.(split(annotation, '_')))))
            end

            hla_types
        end

        data = Escape.HLAData(name, alignment_file, hla_types, tree, missing,
            allele_depth = allele_depth, replacement_mincount = replacement_mincount)
        stan_input = Escape.stan_input(Escape.HLAModel{4}(), data, 
            allele_depth = allele_depth, replacement_mincount = replacement_mincount)

        return HLAData(name, alignment_file, hla_types, tree, stan_input, 
            allele_depth = allele_depth, replacement_mincount = replacement_mincount)
    end
end

function check_inputs(alignment_file, tree_file)
    alignment_error_string = "Alignment file $alignment_file does not exist."
    tree_error_string  = "Phylogenetic tree file $tree_file does not exist."

    isfile(alignment_file) || throw(ArgumentError(alignment_error_string))
    isfile(tree_file) || throw(ArgumentError(tree_error_string))

    r = records(alignment_file)
    tree = phylogenetic_tree(readline(tree_file))
    leaf_names = Set(map(x -> get_property(tree, x, :name), leaves(tree)))
    
    if length(leaf_names) != length(r)
        throw(ArgumentError("Number of leaves of the phylogenetic tree " * 
            "($(length(leaf_names))) does not match number of records " * 
            "($(length(r)))"))
    end

    if Set(leaf_names) != Set(string.(1:length(r))) 
        throw(ArgumentError("Leaves of the phylogenetic tree are not labeled with " *
            "1, 2, ..., $(length(r))."))
    end
    
    for record in r
        annotation = FASTA.identifier(record)
        if length(split(annotation, '_')) != 6
            throw(ArgumentError("Sequence with annotation $annotation does not provide
                6 HLA alleles (use 'missing' to indicate missing alleles)."))
        end
    end

    return nothing
end

name(data::AbstractHLAData) = getfield(data, :name)
records(data::AbstractHLAData) = getfield(data, :records)
hla_types(data::AbstractHLAData) = getfield(data, :hla_types)
stan_input(data::AbstractHLAData) = getfield(data, :stan_input)
allele_depth(data::AbstractHLAData) = getfield(data, :allele_depth)
replacement_mincount(data::AbstractHLAData) = getfield(data, :replacement_mincount)

function compute_stan_input(
    data::AbstractHLAData;
    allele_depth::Int = 1, replacement_mincount::Int = 1
)
    hla_types = Escape.hla_types(data)

    X = Float64.(hla_matrix(hla_types; allele_depth = allele_depth))
    size(X, 2) != 0 || error("Could not parse any HLA alleles with " * 
        "allele_depth = $allele_depth. " *
        "Try setting allele_depth = 1 for HLA alleles with 2-digit accuracy.")
    no_hlas = map(x -> all(x .== 0), eachrow(X))
    no_hlas_idx = findall(no_hlas)
    if any(no_hlas)
        @warn "No HLA alleles parsed for sequences with index: " *
            "$no_hlas_idx. " *
            "These sequences are not included in the model."
    end

    with_hla = .!no_hlas
    r = replacements(data)
    phy = phylogeny_information(data, r)
    Z = epitope_information(data, r, allele_depth)

    y = Int[]
    rs = Int[]
    idx = Int[]

    for i in 1:length(r)
        y_i = targets(r[i], data)
        c = .!ismissing.(y_i) .& with_hla    
        append!(y, y_i[c])
        append!(rs, repeat([i], sum(c)))
        append!(idx, findall(c))
    end

    d = Dict(
        "N" => length(y),
        "S" => size(X, 1),
        "D" => size(X, 2),
        "R" => length(r),
        "X" => X,
        "y" => y,
        "rs" => rs,
        "phy" => phy,
        "Z" => Z,
        "idx" => idx
    )

    return d
end

function epitope_information(
    data::AbstractHLAData, r::Vector{Replacement}, 
    allele_depth::Int
)
    N = length(unique_alleles(data.hla_types, allele_depth = allele_depth))
    Z = Escape.epitope_feature_matrix(
        data, allele_depth = allele_depth
    )[map(x -> x.position, r), :]

    return Z
end

function phylogeny_information(
    data::AbstractHLAData, r::Vector{Replacement}
)
    R = length(r)
    N = length(data.hla_types)
    tree = phylogenetic_tree(data, verbose = true)

    phy = @showprogress "phylogeny " @distributed hcat for i in 1:R
        phylogeny_probabilities(r[i], data, tree)
    end
    phy = phy'

    return phy
end

function is_valid(input::Dict{String, Any})
    m = ErrorException("stan_input is not in a valid format.")

    required_keys = ("Z", "rs", "S","idx", "X", "y", "N", "D", "R", "phy")
    all(map(x -> x in keys(input), required_keys)) || return false
    size(input["Z"]) == (input["R"], input["D"]) || return false
    maximum(input["rs"]) <= input["R"] || return false
    maximum(input["idx"]) <= input["S"] || return false
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