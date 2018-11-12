export HLAAllele, HLAType, parse_allele, is_valid_allele, limit_hla_accuracy, hla_accuracy

struct HLAAllele
    gene::T where T <: AbstractString
    field_1::T where T <: AbstractString
    field_2::Union{Missing, T} where T <: AbstractString
    field_3::Union{Missing, T} where T <: AbstractString
    field_4::Union{Missing, T} where T <: AbstractString
    suffix::Union{Missing, T} where T <: AbstractString
end

struct HLAType
    alleles::NTuple{6, Union{HLAAllele, Missing}}

    function HLAType(alleles::NTuple{6, Union{HLAAllele, Missing}})
        genes = tuple(ismissing(allele) ? missing : allele.gene for allele in alleles)

        for gene in ("A", "B", "C")
            c = count(genes .== gene) 
            if c > 2
                error("HLA type must contain at most 2 $gene alleles, got $(c)")
            end
        end

        new(alleles)
    end
end

function HLAType(::Missing)
    HLAType((missing, missing, missing, missing, missing, missing))
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

        if allele == limit_hla_accuracy(a, depth = hla_accuracy(allele))
            return true
        end
    end

    return false
end

"""
valid HLA alleles must use ':' to seperate fields and have the 'HLA-' prefix.
In contrast to the official nomenclature, 2-digit alleles are allowed 
(e.g. HLA-A*11)
"""
function HLAAllele(s::AbstractString)
    if is_valid_allele(s)
        regex = raw"HLA-([ABC])\*(\d{2}(?!\d)):?((?<=:)\d{2,3}(?!\d))?:?" *
                raw"((?<=:)\d\d(?!\d))?:?((?<=:)\d\d(?!\d))?([NLSCAQ])?" |> Regex

        captures = match(regex, s).captures
        return HLAAllele((i == nothing ? missing : String(i) for i in captures)...)
    else
        throw(error("$s does not follow HLA nomenclature " *
              "(http://hla.alleles.org/nomenclature/naming.html)"))
    end
end

function is_valid_allele(s::AbstractString)
    regex = raw"HLA-([ABC])\*(\d{2}(?!\d)):?((?<=:)\d{2,3}(?!\d))?:?" *
            raw"((?<=:)\d\d(?!\d))?:?((?<=:)\d\d(?!\d))?([NLSCAQ])?" |> Regex

    if match(regex, s) == nothing || length(match(regex, s).captures) < 3
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

function Base.rand(S::Type{HLAAllele}, gene::Symbol)
    gene ∉ [:A, :B, :C] && error("gene must be either :A, :B, or :C")
    fasta_path = joinpath(@__DIR__, "..", "data", "HLA_alleles", 
                         "$(string(gene))_nuc.fasta")
    fasta_reader = BioSequences.FASTA.Reader(open(fasta_path, "r"))
    allele_strings = Vector{String}()

    for record in fasta_reader
        allele_string = split(BioSequences.FASTA.description(record))[1]
        push!(allele_strings, allele_string)
    end

    return parse_allele(rand(allele_strings))
end

function Base.rand(::Type{HLAType})
    hla_alleles = Tuple(rand(HLAAllele, s) for s in (:A, :A, :B, :B, :C, :C))

    return HLAType(hla_alleles)
end

function Base.rand(::Type{HLAType}, d::Int)
    return [rand(HLAType) for i in 1:d]
end

function unique_alleles(hla_types::Vector{HLAType}; depth::Int = 1)
    alleles = Vector{HLAAllele}()
    
    for hla_type in hla_types
        for allele in hla_type.alleles
            hla_type = limit_hla_accuracy(allele, depth = depth)

            if hla_type ∉ alleles
                push!(alleles, hla_type)
            end
        end
    end

    return alleles
end

function limit_hla_accuracy(s::HLAAllele; depth::Int = 1)
    1 <= depth <= 5 || error("depth must be between 1 and 5")

    HLA_components::Vector{Union{Missing, String}} = [missing for i in 1:6]
    HLA_components[1] = s.gene
    fields = fieldnames(HLAAllele)
    for i in 2:(depth + 1)
        HLA_components[i] = getfield(s, fields[i]) 
    end

    return HLAAllele(HLA_components...)
end

function hla_accuracy(s::HLAAllele)
    for (i, field) in enumerate(fieldnames(HLAAllele)[3:end])
        if ismissing(getfield(s, field))
            return i
        end
    end

    return 5
end