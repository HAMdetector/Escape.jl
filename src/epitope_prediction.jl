const NETMHC_ROOTDIR = joinpath(@__DIR__, "..", "data", "netMHC-4.0")

function epitope_feature_matrix(
    data::AbstractHLAData; 
    depth::Int = 1, rank_threshold::Real = 2
)
    df = epitope_prediction(data, rank_threshold = rank_threshold)
    df[!, :allele] = limit_hla_accuracy.(df[!, :allele], depth = depth)
    
    hla_types = Escape.hla_types(data)
    records = Escape.records(data)

    alleles = sort(unique_alleles(hla_types, depth = depth))
    positions = sequence_length(records)

    m = zeros(positions, length(alleles))

    for (i, allele) in enumerate(alleles)
        filtered_df = filter(x -> x[:allele] == allele, df)
        start_positions = filtered_df[!, :start_position]
        stop_positions = filtered_df[!, :stop_position]

        for (start, stop) in zip(start_positions, stop_positions)
            m[start:stop, i] .= 1
        end
    end

    return m
end

function epitope_prediction(data::AbstractHLAData; rank_threshold::Real = 2)
    records = Escape.records(data)
    consensus = consensus_sequence(records)

    df = epitope_prediction(consensus, rank_threshold = rank_threshold)

    return df
end

function epitope_prediction(epitope::String; rank_threshold::Real = 100.0
)
    temp_input = tempname()
    temp_output = tempname()
    # consensus = consensus_sequence(records)

    io = open(temp_input, "w")
    write(io, ">epitope\n")
    write(io, string(replace(epitope, "X" => "-")...))
    close(io) 

    alleles = let
        netmhc_call = read(`$NETMHC_ROOTDIR/netMHC -listMHC -tdir $(tempdir())`, String)
        alleles = split(netmhc_call, '\n')
        filter!(x -> startswith(x, "HLA"), alleles)

        s = ""
        for allele in alleles
            s = s * "$allele,"
        end

        s[1:end-1]
    end

    Base.run(pipeline(`$NETMHC_ROOTDIR/netMHC -tdir $(tempdir()) $temp_input 
        -a $alleles -l 8,9,10,11`, stdout = temp_output))
    df = parse_netmhc(temp_output, rank_threshold = rank_threshold)
    dropmissing!(df)

    rm(temp_output)
    rm(temp_input)

    return sort(df, :rank)
end

function parse_netmhc(filepath::String; rank_threshold::Real = 100.0)
    df = DataFrame(
        [Union{Missing, T} for T in [HLAAllele, String, Int, Int, Float64, Float64]],
        [:allele, :epitope, :start_position, :stop_position, :affinity, :rank]
    )
    for line in readlines(filepath)
        if !any(occursin.(["pos", "#", "Protein", "---", "Linux"], line)) & (line != "")
            components = split(line, r"\s+")

            allele = parse_allele(components[3])
            epitope = components[5]
            start_position = parse(Int, components[2]) + 1
            stop_position = let
                start_position + length(epitope) - length(findall("-", epitope)) - 1
            end
            affinity = parse(Float64, components[14])
            rank = parse(Float64, components[15])
        
            if rank <= rank_threshold
                push!(df, [allele, epitope, start_position, stop_position, affinity, rank])
            end
        end
    end

    return df
end

function consensus_sequence(records::Vector{FASTX.FASTA.Record})
    consensus = Vector{String}(undef, sequence_length(records))

    for i in 1:sequence_length(records)
        sequence = [string(FASTA.sequence(record)[i]) for record in records]
        counts = countmap(sequence)
        consensus[i] = findfirst(x -> x == maximum(values(counts)), counts)
    end

    return string.(consensus...)
end

function sequence_length(records::Vector{FASTX.FASTA.Record})
    sequence_length = maximum([length(FASTA.sequence(r)) for r in records])

    return sequence_length
end