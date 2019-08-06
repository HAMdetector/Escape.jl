const netmhc_rootdir = joinpath(@__DIR__, "..", "data", "netMHC-4.0")

function epitope_prediction(replacement::Replacement, data::HLAData)
    position = replacement.position
    fasta_file = data.fasta_file

    return epitope_prediction(fasta_file, position)
end

function epitope_prediction(filepath::String, position::Int)
    epitope_start = max(1, position - 8)

    epitopes = []
    reader = BioSequences.FASTA.Reader(open(filepath, "r"))
    
    for record in reader
        sequence = BioSequences.FASTA.sequence(record)
        epitope_stop = min(position + 8, length(sequence))

        epitope = string(sequence[epitope_start:epitope_stop])
        epitope = replace(epitope, 'X' => '-')

        push!(epitopes, epitope)
    end

    close(reader)

    consensus = let
        consensus_symbol(x) = sort(collect(countmap(x)), by = x->x[2], rev = true)[1][1]
        epitope_symbols(epitopes, i) = filter(x -> x != '-', [x[i] for x in epitopes])
        c = [consensus_symbol(epitope_symbols(epitopes, i)) for i in 1:length(epitopes[1])]
        string(c...)
    end

    return epitope_prediction(consensus)
end

function epitope_prediction_fasta(filepath::String)
    temp_output = tempname()

    alleles = let
        alleles = split(read(`$netmhc_rootdir/netMHC -listMHC`, String), "\n")
        filter!(x -> startswith(x, "HLA"), alleles)

        s = ""
        for allele in alleles
            s = s * "$allele,"
        end

        s[1:end-1]
    end

    Base.run(pipeline(`$netmhc_rootdir/netMHC $filepath -a $alleles`, stdout = temp_output))
    df = parse_netmhc(temp_output)
    dropmissing!(df)
    
    rm(temp_output)

    return df
end

function epitope_prediction(epitope::String)
    temp_input = tempname()
    temp_output = tempname()

    io = open(temp_input, "w")
    write(io, epitope)
    close(io)

    alleles = let
        alleles = split(read(`$netmhc_rootdir/netMHC -listMHC`, String), "\n")
        filter!(x -> startswith(x, "HLA"), alleles)

        s = ""
        for allele in alleles
            s = s * "$allele,"
        end

        s[1:end-1]
    end

    Base.run(pipeline(`$netmhc_rootdir/netMHC -p $temp_input -a $alleles`, 
                      stdout = temp_output))
    df = parse_netmhc(temp_output)
    dropmissing!(df)

    rm(temp_output)
    rm(temp_input)

    return sort(df, :rank)
end

function parse_netmhc(filepath::String)
    df = DataFrame([Union{Missing, T} for T in [HLAAllele, String, Int, Float64, Float64]],
                   [:allele, :epitope, :offset, :affinity, :rank])
    for line in readlines(filepath)
        if !any(occursin.(["pos", "#", "Protein", "---"], line)) & (line != "")
            components = split(line, r"\s+")
            allele = parse_allele(components[3])
            epitope = components[5]
            offset = parse(Int, components[6])
            affinity = parse(Float64, components[14])
            rank = parse(Float64, components[15])
            
            push!(df, [allele, epitope, offset, affinity, rank])
        end
    end

    return df
end