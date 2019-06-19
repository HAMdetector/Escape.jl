export HLAModel, HLAModelResult, HLAPhylogenyModel, stan_input, classification_accuracy

abstract type HLAModel end
abstract type HLAModelResult end

abstract type HLAPhylogenyModel <: HLAModel end

struct GenericHLAModelResult{T<:HLAModel} <: HLAModelResult
    model::T
    sf::Stanfit
    alleles::Vector{HLAAllele}
    replacement::Replacement
end

function run(model::HLAModel, data::AbstractHLAData, replacement::Replacement;
             iter::Int = 2000, chains::Int = 4, wp::WorkerPool = WorkerPool(), 
             depth::Int = 1)
    
    input = stan_input(model, data, replacement, depth = depth)
    path = model_path(model)

    sf = stan(path, input, iter = iter, chains = chains, wp = wp, 
              stan_args = "adapt delta=0.97")

    # filter relevant parameters to save space
    keep = ["phylogeny_coefficient", "beta_hla", "log_lik", "y_rep"]
    for r in sf.result
        for k in keys(r)
            if !any(startswith.(k, keep))
                delete!(r, k)
            end
        end
    end

    alleles = sort(unique_alleles(filter(x -> missing ∉ x, data.hla_types), depth = depth))

    return GenericHLAModelResult(model, sf, alleles, replacement)
end

function model_path(model::HLAModel)
    error("method model_path(::$(typeof(model))) not defined.")
end

function stan_input(model::HLAPhylogenyModel, data::AbstractHLAData,
                    replacement::Replacement; depth::Int = 1)

    y = targets(replacement, data)
    m = hla_matrix(data.hla_types, depth = depth)

    if ismissing(data.tree)
        ctree = phylogenetic_tree(data)
    else
        ctree = deepcopy(data.tree)
    end

    annotate!(ctree, data, replacement)

    p = state_probabilities(ctree, TwoStateGTR)
    phylogeny_effect = [p[s]["1"] for s in string.(1:length(y))]
    phylogeny_effect = [min(max(0.01, x), 0.99) for x in phylogeny_effect]

    stan_input = Dict("y" => collect(skipmissing(y)), "hla_matrix" => m[.!ismissing.(y), :],
                      "n_entries" => length(collect(skipmissing(y))), 
                      "n_alleles" => size(m)[2], 
                      "phylogeny_effect" => phylogeny_effect[.!ismissing.(y)])
end

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
    filter!(x -> !(0.05 < count(x.second .>= 0) / length(x.second) < 0.95), posteriors)

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