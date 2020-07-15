function HLADataset(s::String)
    HLADataset(Val(Symbol(s)))
end

function HLADataset(::Val{:Rousseau})
    hla_annotation = rousseau_hla_types()

    function extract_hla_types(file::String)
        hla_types = Vector{HLAType}()

        for record in FASTA.Reader(open(file, "r"))
            identifier = FASTA.identifier(record)
            patient_id = split(identifier, '_')[1]
            patient_id = endswith(patient_id, '.') ? patient_id[1:end-1] : patient_id
            accession = length(split(identifier, '_')) == 1 ? "" : split(identifier, '_')[2]
            
            if patient_id in keys(hla_annotation)
                hla_type = hla_annotation[patient_id]
            elseif accession in keys(hla_annotation)
                hla_type = hla_annotation[accession]
            else
                hla_type = HLAType(missing)
            end

            push!(hla_types, hla_type)
        end

        return hla_types
    end

    gag_file = joinpath(dirname(@__DIR__), "data", "test_sequences", "gag_rousseau.fasta")
    gag_hla_types = extract_hla_types(gag_file)
    gag_hla_data = HLAData("gag", gag_file, gag_hla_types, missing, missing)

    nef_file = joinpath(dirname(@__DIR__), "data", "test_sequences", "nef_rousseau.fasta")
    nef_hla_types = extract_hla_types(nef_file)
    nef_hla_data = HLAData("nef", nef_file, nef_hla_types, missing, missing)

    HLADataset("Rousseau", [gag_hla_data, nef_hla_data])
end

function rousseau_hla_types()
    hla_file = joinpath(dirname(@__DIR__), "data", "test_sequences", "hla_rousseau.csv")
    hla_data = readdlm(hla_file, '\t')[2:end, :]

    hla_types = Dict{String, HLAType}()

    for i in 1:size(hla_data, 1)
        patient_id, accession = string.(hla_data[i, 1:2])
        
        hla_alleles = parse_allele.(split(hla_data[i, 3], ' '))
        while length(hla_alleles) < 6 
            push!(hla_alleles, missing) 
        end
        
        hla_type = HLAType(tuple(hla_alleles...))
        if patient_id != "" hla_types[patient_id] = hla_type end
        if accession != "" hla_types[accession] = hla_type end
    end

    return hla_types
end