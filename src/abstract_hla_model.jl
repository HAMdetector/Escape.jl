export HLAModel, HLAModelResult, HLAPhylogenyModel, stan_input, classification_accuracy

abstract type HLAModel end
abstract type HLAModelResult end

abstract type HLAPhylogenyModel <: HLAModel end


function stan_input(model::HLAModel, data::AbstractHLAData, replacement::Replacement;
                    depth::Int = 1)
    y = targets(replacement, data)
    m = hla_matrix(data.hla_types, depth = depth)

    stan_input = Dict("y" => collect(skipmissing(y)), "hla_matrix" => m[.!ismissing.(y), :],
                      "n_entries" => length(collect(skipmissing(y))), 
                      "n_alleles" => size(m)[2])
end

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
    
    d = Dict{HLAAllele, Vector{Float64}}()
    for (i, allele) in enumerate(result.alleles)
        d[allele] = posterior["beta_hla.$i"]
    end

    return d
end

function posterior_interval(v::Vector{T}; cutoff::Real = 0.95) where T <: Real
    sorted = sort(v)
    i1 = round(Int, (1 - cutoff) * length(v))
    i2 = round(Int, cutoff * length(v))

    return (sorted[i1], sorted[i2])
end