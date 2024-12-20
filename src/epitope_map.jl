export EpitopeMap

struct EpitopeMap 
    name::String
    epitopes::Vector{String}
    starts::Vector{Int}
    stops::Vector{Int}
    alleles::Vector{HLAAllele}
    maplength::Int

    function EpitopeMap(name, epitopes, starts, stops, alleles, maplength)
        if !(length(epitopes) == length(starts) == length(stops))
            error("epitope sequence and start and stop positions" * 
                  "must have the same length.")
        elseif !all(starts .<= stops)
            error("epitope start must be less than epitope stop.")
        elseif !(all((starts .> 0) .& (stops .> 0)))
            error("start and stop positions must be positive.")
        elseif !all(length.(epitopes) .== (stops .- starts .+ 1))
            error("epitope length must match the given start and stop positions.")
        elseif (length(stops) != 0) && !(maplength >= maximum(stops))
            error("maplength must be >= largest stop position.") 
        end

        new(name , epitopes, starts, stops, alleles, maplength)
    end
end

function epitope_map(data::AbstractHLAData; rank_threshold::Real = 100)
    df = epitope_prediction(data, rank_threshold = rank_threshold)
    
    epitopes = df[!, :peptide]
    starts = df[!, :pos] .+ 1
    stops = starts .+ length.(df[!, :peptide]) .- 1
    alleles = parse_allele.(df[!, :best_allele])
    maplength = sequence_length(data)

    map = EpitopeMap(name(data), epitopes, starts, stops, alleles, maplength)
    
    return map
end

function Base.in(x::Pair{Replacement, HLAAllele}, y::EpitopeMap)
    for (i, (start, stop)) in enumerate(zip(y.starts, y.stops))
        if (start <= x[1].position <= stop) & (y.alleles[i] == x[2])
            return true
        end
    end

    return false
end

function match_epitope(epitope::String, data::AbstractHLAData; 
        max_levenshtein_distance::Int = 0
)
    position = Int[]
    query = ApproximateSearchQuery(LongAA(epitope))

    for record in records(data)
        sequence = FASTX.sequence(record)
        sequence_string = replace(string(sequence), 'X' => '-')
        sequence = LongAA(sequence_string)

        idx = findfirst(query, max_levenshtein_distance, sequence)

        if idx != nothing
            push!(position, idx[1])
        end
    end

    length(unique(position)) == 1 || return nothing

    return position[1]
end
