const NETMHC_ROOTDIR = ENV["NETMHC_ROOTDIR"]

function epitope_feature_matrix(
    data::AbstractHLAData; 
    depth::Int = 1, rank_threshold::Real = 2,
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
        start_positions = filtered_df[!, :position]
        stop_positions = filtered_df[!, :position] .+ length.(filtered_df[!, :peptide]) .- 1

        for (start, stop) in zip(start_positions, stop_positions)
            m[start:stop, i] .= 1
        end
    end

    return m
end

function epitope_prediction(data::AbstractHLAData; rank_threshold::Real = 2)
    records = Escape.records(data)
    consensus = consensus_sequence(records)

    alleles = unique_alleles(hla_types(data))

    df = epitope_prediction(consensus, alleles, rank_threshold = rank_threshold)

    return df
end

function epitope_prediction(
    query::String, alleles::Vector{HLAAllele};
    rank_threshold::Real = 100
)
    netmhc_alleles = Escape.netmhc_alleles()
    relevant_alleles = HLAAllele[]
    
    
    for allele in alleles
        if hla_accuracy(allele) == 1
            for netmhc_allele in netmhc_alleles
                if limit_hla_accuracy(netmhc_allele, depth = 1) == allele
                    push!(relevant_alleles, netmhc_allele)
                end
            end
        else
            push!(relevant_alleles, limit_hla_accuracy(allele, depth = 2))
        end
    end

    relevant_alleles = sort(unique(relevant_alleles))
    allele_partitions = collect.(Base.Iterators.partition(relevant_alleles, 10))

    p = Progress(length(allele_partitions), desc = "epitope prediction ")

    dfs = asyncmap(allele_partitions, ntasks = Threads.nthreads()) do x
        next!(p)
        epitope_prediction_limited(query, x, rank_threshold = rank_threshold)
    end

    combined_df = dfs[1]
    for i in 2:length(dfs)
        append!(combined_df, dfs[i])
    end

    filter!(x -> x[:rank_percent] <= rank_threshold, combined_df)
    sort!(combined_df, :rank_percent)

    return combined_df
end

function epitope_prediction_limited(
    query::String, alleles::Vector{HLAAllele};
    rank_threshold::Real = 100
)
    if (length(alleles) > 85)
        error("maxmimum number of 85 alleles supported, got $(length(alleles)).")
    elseif any(ismissing.(getfield.(alleles, :field_2)))
        error("only ")
    end

    length(alleles) <= 85 & all(.!ismissing.(getfield.(alleles, :field_2)))

    temp_input = tempname()
    temp_output = tempname()

    open(temp_input, "w") do io
        write(io, ">epitope\n")
        write(io, string(replace(query, "X" => "-")...))
    end

    alleles = sort(unique(limit_hla_accuracy.(alleles, depth = 2)))

    s = ""
    for allele in alleles
        if ismissing(allele.field_2)
            error("field 2 is missing.")
        else
            s = s * "HLA-" * allele.gene * allele.field_1 * allele.field_2 * ','
        end
    end
    s = s[1:end-1]

    mhc_output = Base.read(`$NETMHC_ROOTDIR/netMHCpan $temp_input -a $s`) |> String
    df = parse_netmhc(mhc_output)
    filter!(x -> x[:rank_percent] <= rank_threshold, df)

    return df
end

function parse_netmhc(output::String)
    df = DataFrame(
        position = Int[], allele = HLAAllele[], peptide = String[], core = String[],
        core_start = Int[], deletion_position = Int[], deletion_length = Int[],
        insertion_position = Int[], insertion_length = Int[],
        interaction_core = String[], identity = String[], score = Float64[],
        rank_percent = Float64[]
    )

    for line in split(output, '\n')
        if occursin("HLA-", line) & !(any(startswith.(line, ['#', "HLA-", "Protein"])))
            components = split(line, r"\s+")
            components = string.(components[2:end])

            position = parse(Int, components[1])
            allele = parse_allele(components[2])
            peptide = components[3]
            core = components[4]
            core_start = parse(Int, components[5])
            deletion_position = parse(Int, components[6])
            deletion_length = parse(Int, components[7])
            insertion_position = parse(Int, components[8])
            insertion_length = parse(Int, components[9])
            interaction_core = components[10]
            identity = components[11]
            score = parse(Float64, components[12])
            rank_percent = parse(Float64, components[13])

            push!(df, [position, allele, peptide, core, core_start, deletion_position,
                deletion_length, insertion_position, insertion_length,
                interaction_core, identity, score, rank_percent])
        end
    end

    return df
end

function netmhc_alleles()
    netmhc_call = read(`$NETMHC_ROOTDIR/netMHCpan -listMHC`, String)
    alleles = split(netmhc_call, '\n')

    filter!(x -> startswith(x, "HLA"), alleles)
    filter!(x -> !startswith(x, "HLA-E"), alleles)
    filter!(x -> !startswith(x, "HLA-G"), alleles)
    filter!(x -> !(':' in x), alleles)

    alleles = parse_allele.(alleles)

    return alleles
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