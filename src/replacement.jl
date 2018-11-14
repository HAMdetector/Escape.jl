export Replacement, replacements

struct Replacement
    protein::String
    position::Int
    replacement::Char
end

function replacements(data::HLAData; mincount::Int = 5)
    replacements = Vector{Replacement}()

    reader = BioSequences.FASTA.Reader(open(data.fasta_file, "r"))
    N = minimum(length(BioSequences.FASTA.sequence(record)) for record in reader)

    for i in 1:N
        reader = BioSequences.FASTA.Reader(open(data.fasta_file, "r"))
        counts = [Char(BioSequences.FASTA.sequence(record)[i]) for record in reader] |>
            StatsBase.countmap 
        filter!(x -> x.first ∉ ('-', 'X') && x.second >= mincount, counts)
        
        if length(counts) > 1
            for (k, v) in counts
                push!(replacements, Replacement(data.name, i, k))
            end
        end
    end

    close(reader)
    
    return replacements
end