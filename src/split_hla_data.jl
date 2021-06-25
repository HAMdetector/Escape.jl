struct SplitHLAData <: AbstractHLAData
    name::String
    records::Vector{FASTX.FASTA.Record}
    hla_types::Vector{HLAType}
    tree::Union{Missing, PhylogeneticTree}
    stan_input::Dict{String, Any}
    original_data::HLAData
    split::UnitRange{Int64}
    split_idx::Int

    function SplitHLAData(
        name::String, 
        records::Vector{FASTX.FASTA.Record}, 
        hla_types::Vector{HLAType},
        tree::Union{Missing, PhylogeneticTree},
        stan_input::Union{Missing, Dict{String, Any}},
        original_data::HLAData,
        split::UnitRange{Int64},
        split_idx::Int
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
        
        new(name, records, hla_types, tree, stan_input, original_data, split, split_idx)
    end
end

function split_hla_data(data::HLAData, n::Int = 10)
    split_data = SplitHLAData[]
    replacements = Escape.replacements(data, mincount = 1)
    stan_input = Escape.stan_input(HLAModel{4}(), data)
    R = stan_input["R"]
    splits = Iterators.partition(1:R, ceil(Int, R // n))

    for (i, split) in collect(enumerate(splits))
        stan_input_split = copy(stan_input)
        in_split = findall(x -> x in split, stan_input_split["rs"])

        stan_input_split["Z"] = stan_input_split["Z"][split, :]
        stan_input_split["rs"] = map(x -> findfirst(y -> y == x, split), 
            stan_input_split["rs"][in_split])
        stan_input_split["idx"] = stan_input_split["idx"][in_split]
        stan_input_split["y"] = stan_input_split["y"][in_split]
        stan_input_split["N"] = length(in_split)
        stan_input_split["R"] = length(split)
        stan_input_split["phy"] = stan_input_split["phy"][split, :]
        
        split_data_ = SplitHLAData(data.name, data.records, data.hla_types, data.tree,
            stan_input_split, data, split, i)
        
        @assert length(Escape.replacements(split_data_)) == stan_input_split["R"]
        
        push!(split_data, split_data_)
    end

    return split_data
end