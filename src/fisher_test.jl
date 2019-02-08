export FisherTest, FisherTestResult

struct FisherTest end

struct FisherTestResult 
    counts::Dict{HLAAllele, Array{Int, 2}}
    p_values::Dict{HLAAllele, Float64}
    log_odds::Dict{HLAAllele, Float64}
end

function Escape.run(model::FisherTest, data::AbstractHLAData, replacement::Replacement)
    y = targets(replacement, data)
    m = hla_matrix(data.hla_types)
    alleles = sort(unique_alleles(filter(x -> missing ∉ x, data.hla_types)))

    counts = Dict{HLAAllele, Array{Int, 2}}()
    p_values = Dict{HLAAllele, Float64}()
    log_odds = Dict{HLAAllele, Float64}()

    for (i, allele) in enumerate(alleles)
        mutation_present = Bool.(collect(skipmissing(y)))
        allele_present = Bool.(min.(1, m[.!ismissing.(y), i]))

        c, p, beta = fisher_exact_test(mutation_present, allele_present)
        counts[allele] = c
        p_values[allele] = p
        log_odds[allele] = beta
    end

    return FisherTestResult(counts, p_values, log_odds)
end

function fisher_exact_test(a::AbstractVector{Bool}, b::AbstractVector{Bool})
    c = zeros(Int, 2, 2)

    for i in eachindex(a)
        c[2 - a[i], 2 - b[i]] += 1
    end

    if (c[:, 1] == [0, 0]) | (c[:, 2] == [0, 0])
        return (counts = c, p = 1, beta = Inf)
    end

    f = FisherExactTest(c[1,1], c[1,2], c[2,1], c[2,2])

    return (counts = c, p = pvalue(f), beta = log(f.ω))
end