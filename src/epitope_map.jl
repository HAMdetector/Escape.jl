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
        elseif !all(length.(epitopes) .- count.(x -> x == '-', epitopes) .== 
                (stops .- starts .+ 1))
            error("epitope length must match the given start and stop positions.")
        elseif !(maplength >= maximum(stops))
            error("maplength must be >= largest stop position.") 
        end

        new(name , epitopes, starts, stops, alleles, maplength)
    end
end

function epitope_map(data::AbstractHLAData; rank_threshold::Real = 100)
    df = epitope_prediction_fasta(data.fasta_file, rank_threshold = rank_threshold)
    
    epitopes = df[!, :epitope]
    starts = df[!, :start_position]
    stops = df[!, :stop_position]
    alleles = df[!, :allele]
    maplength = fasta_length(data.fasta_file)

    map = EpitopeMap(data.name, epitopes, starts, stops, alleles, maplength)
    
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