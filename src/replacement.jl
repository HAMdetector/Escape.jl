export Replacement, replacements

struct Replacement
    protein::String
    position::Int
    replacement::Char
    negated::Bool
end

protein(r::Replacement) = getfield(r, :protein)
position(r::Replacement) = getfield(r, :position)
replacement(r::Replacement) = getfield(r, :replacement)
negated(r::Replacement) = getfield(r, :negated)

function replacement(
        protein::String, position::Int, replacement::Char; 
        negated::Bool = false
    )

    return Replacement(protein, position, replacement, negated)
end

function replacements(data::AbstractHLAData; mincount::Int = 1)
    replacements = Replacement[]
    records = Escape.records(data)
    N = minimum(length(FASTA.sequence(record)) for record in records)

    for i in 1:N
        counts = [Char(FASTA.sequence(record)[i]) for record in records] |>
            StatsBase.countmap 
        filter!(x -> x.first ∉ ('X', '*', '-') && x.second >= mincount, counts)
        
        if length(counts) > 1
            for (k, v) in counts
                push!(replacements, replacement(name(data), i, k))
            end
        end
    end
    
    return replacements
end

function targets(replacement::Replacement, data::AbstractHLAData)
    t = Vector{Union{Missing, Int}}()

    for record in records(data)
        r = Escape.replacement(replacement)
        p = position(replacement)

        symbol = Char(FASTA.sequence(record)[p])
        r = replacement.replacement
        
        if negated(replacement)
            if symbol == r
                push!(t, 0)
            elseif symbol in ('X', '-', '*')
                push!(t, missing)
            else
                push!(t, 1)
            end
        elseif !negated(replacement)
            if symbol == r
                push!(t, 1)
            elseif symbol in ('X', '-', '*')
                push!(t, missing)
            else
                push!(t, 0)
            end
        end
    end

    return t
end