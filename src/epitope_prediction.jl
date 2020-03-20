const NETMHC_ROOTDIR = joinpath(@__DIR__, "..", "data", "netMHC-4.0")

function epitope_feature_matrix(
    data::AbstractHLAData; 
    depth::Int = 1, rank_threshold::Real = 2
)
    df = epitope_prediction(data, rank_threshold = rank_threshold)
    df[!, :allele] = limit_hla_accuracy.(df[!, :allele], depth = depth)
    
    alleles = sort(unique_alleles(data.hla_types, depth = depth))
    positions = fasta_length(data.fasta_file)

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
    df = epitope_prediction_fasta(data.fasta_file, rank_threshold = rank_threshold)

    return df
end

function epitope_prediction(epitope::String; rank_threshold::Real = 100.0)
    temp_input = tempname()
    io = open(temp_input, "w")
    write(io, ">epitope\n")
    write(io, replace(epitope, "X" => "-"))
    close(io)

    df = epitope_prediction_fasta(temp_input, rank_threshold = rank_threshold)
    rm(temp_input)

    return df
end

function epitope_prediction_fasta(filepath::String; rank_threshold::Real = 100.0)
    temp_input = tempname()
    temp_output = tempname()
    consensus = consensus_sequence(filepath)

    io = open(temp_input, "w")
    write(io, ">consensus\n")
    write(io, string(replace(consensus, "X" => "-")...))
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

function consensus_sequence(filepath::String)
    consensus = Vector{String}(undef, fasta_length(filepath))
    reader = FASTA.Reader(open(filepath, "r"))

    for i in 1:fasta_length(filepath)
        reader = FASTA.Reader(open(filepath, "r"))
        sequence = [string(FASTA.sequence(record)[i]) for record in reader]
        close(reader)

        counts = countmap(sequence)
        consensus[i] = findfirst(x -> x == maximum(values(counts)), counts)
    end

    close(reader)
    return consensus
end

function fasta_length(filepath::String)
    reader = FASTA.Reader(open(filepath, "r"))
    fasta_length = maximum([length(FASTA.sequence(r)) for r in reader])

    close(reader)

    return fasta_length
end