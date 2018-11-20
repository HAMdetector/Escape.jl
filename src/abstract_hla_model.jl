export HLAModel, HLAModelResult, HLAPhylogenyModel, classification_accuracy

abstract type HLAModel end
abstract type HLAModelResult end

abstract type HLAPhylogenyModel <: HLAModel end

function classification_accuracy(result::HLAModelResult)
    n_entries = result.sf.data["n_entries"]
    posterior = StanInterface.extract(result.sf)
    
    estimated = [median(posterior["y_rep.$i"]) > 0.5 ? 1 : 0 for i in 1:n_entries]
    observed = result.sf.data["y"]
    
    return Float64(count(estimated .== observed) / n_entries)
end

function relevant_alleles(result::HLAModelResult)
    posteriors = allele_posteriors(result)
    filter!(x -> !(posterior_interval(x.second)[1] < 0 <= posterior_interval(x.second)[2]),
            posteriors)

    return posteriors
end

function allele_posteriors(result::HLAModelResult)
    posterior = extract(result.sf)
    
    alleles = Vector{Pair{HLAAllele, Vector{Float64}}}()
    for (i, allele) in enumerate(result.alleles)
        push!(alleles, allele => posterior["beta_hla.$i"])
    end

    return alleles
end

function posterior_interval(v::Vector{T}; cutoff::Real = 0.95) where T <: Real
    sorted = sort(v)
    i1 = round(Int, (1 - cutoff) * length(v))
    i2 = round(Int, cutoff * length(v))

    return (sorted[i1], sorted[i2])
end