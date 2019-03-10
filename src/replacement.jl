export Replacement, replacements

struct Replacement
    protein::String
    position::Int
    replacement::Char
end

function replacements(data::AbstractHLAData; mincount::Int = 10)
    replacements = Vector{Replacement}()

    reader = BioSequences.FASTA.Reader(open(data.fasta_file, "r"))
    N = minimum(length(BioSequences.FASTA.sequence(record)) for record in reader)

    for i in 1:N
        reader = BioSequences.FASTA.Reader(open(data.fasta_file, "r"))
        counts = [Char(BioSequences.FASTA.sequence(record)[i]) for record in reader] |>
            StatsBase.countmap 
        filter!(x -> x.first ∉ ('X', '*') && x.second >= mincount, counts)
        
        if length(counts) > 1
            for (k, v) in counts
                push!(replacements, Replacement(data.name, i, k))
            end
        end
    end

    close(reader)
    
    return replacements
end

function targets(replacement::Replacement, data::AbstractHLAData)
    reader = BioSequences.FASTA.Reader(open(data.fasta_file, "r"))
    t = Vector{Union{Missing, Int}}()

    for record in reader
        symbol = Char(BioSequences.FASTA.sequence(record)[replacement.position])

        if symbol == replacement.replacement
            push!(t, 1)
        elseif symbol ∈ ('X', '-', '*')
            push!(t, missing)
        else
            push!(t, 0)
        end
    end

    close(reader)
    return t
end