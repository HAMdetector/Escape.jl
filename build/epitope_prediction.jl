function epitope_feature_matrix(
    data::AbstractHLAData; rank_threshold::Real = 2, allele_depth::Int = 1,
)
    df = epitope_prediction(data, rank_threshold = rank_threshold, 
        allele_depth = allele_depth)
    
    hla_types = Escape.hla_types(data)
    alleles = sort(unique_alleles(hla_types, allele_depth = allele_depth))
    
    positions = sequence_length(records(data))

    m = zeros(positions, length(alleles))

    for row in eachrow(df)
        allele_idx = findfirst(x -> x == parse_allele(row[:best_allele]), alleles)
        start_idx = row[:pos] + 1
        stop_idx = start_idx + length(row[:peptide]) - 1

        m[start_idx:stop_idx, allele_idx] .= 1
    end

    return m
end

function epitope_prediction(data::AbstractHLAData; 
    rank_threshold::Real = 2, allele_depth::Int = 1
)
    
    records = Escape.records(data)
    consensus = replace(consensus_sequence(records), '-' => 'X')
    alleles = unique_alleles(hla_types(data), allele_depth = allele_depth)

    df = epitope_prediction(consensus, alleles, rank_threshold = rank_threshold)

    return df
end

function epitope_prediction(
    query::String, alleles::Vector{HLAAllele};
    rank_threshold::Real = 100
)
    @assert all(hla_accuracy.(alleles) .== hla_accuracy(alleles[1]))
    accuracy = min(2, hla_accuracy(alleles[1]))

    valid_4_digits = valid_alleles()
    filter!(x -> parse(Int, x.field_2) <= 1, valid_4_digits)

    valid_matching_digits = limit_hla_accuracy.(
        valid_4_digits, 
        allele_depth = accuracy
    )

    alleles = valid_4_digits[map(x -> x in alleles, valid_matching_digits)]
    allele_strings = String[]
    output_file = tempname()

    for allele in string.(alleles)
        push!(allele_strings, allele[1:5] * allele[7:end])
    end
    
    @suppress run(
        `$(Conda.BINDIR)/mhcflurry-predict-scan
            --alleles $allele_strings 
            --sequences $query
            --out $output_file
            --results-all`
    )

    df = DataFrame(CSV.File(output_file))
    filter!(x -> (x[:affinity_percentile] <= rank_threshold) | 
        (x[:presentation_percentile] <= rank_threshold), df)

    df[!, :best_allele] = Escape.limit_hla_accuracy.(
        parse_allele.(df[!, :best_allele]), 
        allele_depth = accuracy
    ) .|> string

    return df
end

function valid_alleles()
    mhcflurry_output = @suppress read(
        `$(Conda.BINDIR)/mhcflurry-predict --list-supported-alleles`, String
    )
    valid_allele_strings = split(mhcflurry_output, '\n')
    filter!(x -> startswith(x, "HLA"), valid_allele_strings)

    valid_alleles = parse_allele.(valid_allele_strings)

    return unique(skipmissing(valid_alleles))
end

function consensus_sequence(records::Vector{FASTX.FASTA.Record})
    consensus = Vector{String}(undef, sequence_length(records))

    for i in 1:sequence_length(records)
        sequence = [string(FASTA.sequence(record)[i]) for record in records]
        counts = countmap(sequence)
        delete!(counts, "-")
        delete!(counts, "X")

        if isempty(counts)
            consensus[i] = "-"
        else
            consensus[i] = findfirst(x -> x == maximum(values(counts)), counts)
        end 
    end

    return string.(consensus...)
end

function sequence_length(records::Vector{FASTX.FASTA.Record})
    sequence_length = maximum([length(FASTA.sequence(r)) for r in records])

    return sequence_length
end