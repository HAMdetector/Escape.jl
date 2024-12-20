export HLAAllele, HLAType, parse_allele, is_valid_allele, limit_hla_accuracy, hla_accuracy

struct HLAAllele
    gene::String
    field_1::String
    field_2::Union{Missing, String}
    field_3::Union{Missing, String}
    field_4::Union{Missing, String}
    suffix::Union{Missing, String}
end

struct HLAType
    alleles::NTuple{6, Union{HLAAllele, Missing}}

    function HLAType(alleles::NTuple{6, Union{HLAAllele, Missing}})
        genes = tuple(ismissing(allele) ? missing : allele.gene for allele in alleles)

        for gene in ("A", "B", "C")
            c = count(genes .== gene) 
            if c > 2
                error("HLA type must contain at  2 $gene alleles, got $(c)")
            end
        end

        new(alleles)
    end
end

function HLAType(::Missing)
    HLAType((missing, missing, missing, missing, missing, missing))
end

function hla_matrix(data::Vector{HLAType}; allele_depth::Int = 1)
    alleles = sort(unique_alleles(data, allele_depth = allele_depth))
    m = zeros(Int, length(data), length(alleles))

    for i in eachindex(data)
        hla_type = data[i]
        for allele in hla_type.alleles
            a = limit_hla_accuracy(allele, allele_depth = allele_depth)
            (ismissing(a) || hla_accuracy(a) != allele_depth) && continue

            m[i, findfirst(x -> x == a, alleles)] = 1
        end
    end

    return m
end

function hla_matrix_standardized(data::Vector{HLAType}; allele_depth::Int = 1)
    alleles = sort(unique_alleles(data, allele_depth = allele_depth))
    m = zeros(Float64, length(data), length(alleles))

    for (i, allele) in enumerate(alleles)
        frequency = count(x -> allele in x, data) / length(data)
        for (j, hla_type) in enumerate(data)
            if allele in hla_type
                m[j, i] = 1 - frequency
            else
                m[j, i] = -frequency
            end
        end
    end

    return m
end

function Base.in(::Missing, hla_type::HLAType)
    if any(ismissing.(hla_type.alleles))
        return true
    end
    
    return false
end

function Base.in(allele::HLAAllele, hla_type::HLAType)
    for a in hla_type.alleles
        ismissing(a) && continue

        if allele == limit_hla_accuracy(a, allele_depth = hla_accuracy(allele))
            return true
        end
    end

    return false
end

function HLAAllele(s::AbstractString)
    if is_valid_allele(s)
        regex = raw"HLA-([ABC])\*(\d{2}(?!\d)):?((?<=:)\d{2,3}(?!\d))?:?" *
                raw"((?<=:)\d\d(?!\d))?:?((?<=:)\d\d(?!\d))?([NLSCAQ])?" |> Regex

        captures = match(regex, s).captures
        return HLAAllele((i === nothing ? missing : String(i) for i in captures)...)
    else
        throw(error("$s does not follow HLA nomenclature " *
              "(http://hla.alleles.org/nomenclature/naming.html)"))
    end
end

function is_valid_allele(s::AbstractString)
    regex = raw"HLA-([ABC])\*(\d{2}(?!\d)):?((?<=:)\d{2,3}(?!\d))?:?" *
            raw"((?<=:)\d\d(?!\d))?:?((?<=:)\d\d(?!\d))?([NLSCAQ])?" |> Regex

    if match(regex, s) === nothing || length(match(regex, s).captures) < 3
        return false
    end

    return true
end

function Base.isless(x::HLAAllele, y::HLAAllele)
    for name in fieldnames(typeof(x))
        x_value = getfield(x, name)
        y_value = getfield(y, name)

        xor(ismissing(x_value), ismissing(y_value)) && return missing
        ismissing(x_value) && ismissing(y_value) && return false

        if x_value < y_value
            return true
        elseif x_value > y_value
            return false
        end
    end

    return false
end

function parse_allele(::Missing)
    return missing
end

function parse_allele(x::AbstractString)
    if occursin("non", x)
        return HLAAllele(missing)
    end

    s = replace(x, r"[()]" => s"")

    # place * after gene, remove w in Cw, add 'HLA-' prefix
    s = replace(s, r"([ABC])w?\*?(\d.*)" => s"HLA-\1*\2")

    # append 0 for single-digit allele groups
    s = replace(s, r"HLA-([ABC])w?\*?([1-9])(\D|$)(.*)" => s"HLA-\1*0\2")

    # place field seperator after allele group
    s = replace(s, r"HLA-([ABC]\*)(\d{2}):?(.+)" => s"HLA-\1\2:\3")

    # place field seperator for 3 digit field_1
    s = replace(s, r"HLA-([ABC]\*\d{2}:)(\d{3}):?(\d{2}):?(\d{2})(.*)" => 
                s"HLA-\1\2:\3:\4\5")

    # place field seperator for 2 digit field_1
    s = replace(s, r"HLA-([ABC]\*\d{2}:)(\d{2}):?(\d{2}):?(\d{2})(.*)" =>
                s"HLA-\1\2:\3:\4\5")

    # place field seperator for 6 digit HLA alleles
    s = replace(s , r"HLA-([ABC]\*\d{2}:)(\d{2}):?(\d{2})(.*)" => s"HLA-\1\2:\3\4")

    if is_valid_allele(s)
        return HLAAllele(s)
    else
        return missing
    end
end

function parse_allele(x...)
    return NTuple{length(x), HLAAllele}(parse_allele(allele) for allele in x)
end

function convert(::Type{String}, x::HLAAllele)
    s = "HLA-"

    s = s * x.gene * "*" * x.field_1
    
    s = !ismissing(x.field_2) ? s * ":" * x.field_2 : s
    s = !ismissing(x.field_3) ? s * ":" * x.field_3 : s
    s = !ismissing(x.field_4) ? s * ":" * x.field_4 : s
    s = !ismissing(x.suffix) ? s * x.suffix : s
    
    return s
end

function Base.show(io::IO, x::HLAAllele)
    print(io, convert(String, x))
end


Base.rand(S::Type{HLAAllele}, gene::Symbol) = Base.rand(Random.default_rng(), S, gene)

function Base.rand(
        rng::AbstractRNG, 
        S::Type{HLAAllele}, 
        gene::Symbol
)
    gene ∉ [:A, :B, :C] && error("gene must be either :A, :B, or :C")
    fasta_path = joinpath(@__DIR__, "..", "data", "HLA_alleles", 
                         "$(string(gene))_nuc.fasta")
    fasta_reader = FASTA.Reader(open(fasta_path, "r"))
    allele_strings = Vector{String}()

    for record in fasta_reader
        allele_string = split(FASTA.description(record), " ")[2]
        push!(allele_strings, allele_string)
    end

    return parse_allele(rand(rng, allele_strings))
end

Base.rand(T::Type{HLAType}) = Base.rand(Random.default_rng(), T)

function Base.rand(rng::AbstractRNG, ::Type{HLAType})
    hla_alleles = Tuple(rand(rng, HLAAllele, s) for s in (:A, :A, :B, :B, :C, :C))

    return HLAType(hla_alleles)
end

Base.rand(T::Type{HLAType}, d::Int) = Base.rand(Random.default_rng(), T::Type{HLAType}, d::Int)

function Base.rand(rng::AbstractRNG, T::Type{HLAType}, d::Int)
    return [rand(rng, HLAType) for i in 1:d]
end

function unique_alleles(hla_types::Vector{HLAType}; allele_depth::Int = 1)
    alleles = Vector{HLAAllele}()
    
    for hla_type in hla_types
        for allele in hla_type.alleles
            a = limit_hla_accuracy(allele, allele_depth = allele_depth)
            (ismissing(a) || hla_accuracy(a) != allele_depth) && continue

            if a ∉ alleles
                push!(alleles, a)
            end
        end
    end

    return sort(alleles)
end

function limit_hla_accuracy(::Missing; allele_depth::Int = 1)
    return missing
end

function limit_hla_accuracy(s::HLAAllele; allele_depth::Int = 1)
    1 <= allele_depth <= 5 || error("allele_depth must be between 1 and 5.")

    HLA_components = Vector{Union{Missing, String}}(undef, 6)
    fill!(HLA_components, missing)

    HLA_components[1] = s.gene
    HLA_components[2] = s.field_1
    
    if allele_depth >= 2
        HLA_components[3] = s.field_2
    elseif allele_depth >= 3
        HLA_components[4] = s.field_3
    elseif allele_depth >= 4
        HLA_components[5] = s.field_4
    elseif allele_depth == 6
        HLA_components[6] = s.suffix
    end

    return HLAAllele(HLA_components[1], HLA_components[2], HLA_components[3],
        HLA_components[4], HLA_components[5], HLA_components[6])
end

function hla_accuracy(s::HLAAllele)
    for (i, field) in enumerate(fieldnames(HLAAllele)[3:end])
        if ismissing(getfield(s, field))
            return i
        end
    end

    return 5
end
