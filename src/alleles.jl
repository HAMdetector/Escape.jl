export HLAAllele, parse_allele, is_valid_allele

struct HLAAllele
    gene::Union{T} where T <: AbstractString
    field_1::Union{T} where T <: AbstractString
    field_2::Union{Missing, T} where T <: AbstractString
    field_3::Union{Missing, T} where T <: AbstractString
    field_4::Union{Missing, T} where T <: AbstractString
    suffix::Union{Missing, T} where T <: AbstractString
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

function HLAAllele(::Missing)
    HLAAllele((missing for i in 1:6)...)
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
    return collect(parse_allele.(x))
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