export BernoulliPhylogenyModel, BernoulliPhylogenyResult

struct BernoulliPhylogenyModel <: HLAPhylogenyModel
    prior::Symbol
    iter::Int
    chains::Int

    function BernoulliPhylogenyModel(prior, iter, chains)
        if prior ∉ (:broad_t, :finnish_horseshoe, :narrow_t)
            error("prior must be one :broad_t or :finnish_horseshoe, got $prior.")
        end

        new(prior, iter, chains)
    end
end

struct BernoulliPhylogenyResult <: HLAModelResult
    sf::Stanfit
    alleles::Vector{HLAAllele}
    replacement::Replacement
end

function BernoulliPhylogenyModel(iter::Int, chains::Int)
    return BernoulliPhylogenyModel(:finnish_horseshoe, iter, chains)
end

function BernoulliPhylogenyModel(; prior = :finnish_horseshoe, iter = 2000, chains = 4)
    return BernoulliPhylogenyModel(prior, iter, chains)
end

function run(model::BernoulliPhylogenyModel, data::AbstractHLAData, 
             replacement::Replacement; wp::WorkerPool = WorkerPool())
    if model.prior == :broad_t
        path = joinpath(@__DIR__, "..", "data", "stan", "bernoulli_phylogeny")
    elseif model.prior == :finnish_horseshoe
        path = joinpath(@__DIR__, "..", "data", "stan", "bernoulli_phylogeny_hs")
    end

    if ismissing(data.tree)
        input = stan_input(model, data, replacement, phylogenetic_tree(data))
    else
        input = stan_input(model, data, replacement, data.tree)
    end
    
    sf = stan(path, input, chains = model.chains, iter = model.iter, wp = wp,
              stan_args = "adapt delta=0.97")
    
    # filter relevant parameters to save space
    keep = ["intercept", "phylogeny_coefficient", "beta_hla", "log_lik", "y_rep",
            "theta", "lp__"]
    for r in sf.result
        for k in keys(r)
            if !any(startswith.(k, keep))
                delete!(r, k)
            end
        end
    end

    alleles = sort(unique_alleles(data.hla_types))

    return BernoulliPhylogenyResult(sf, alleles, replacement)
end


function stan_input(model::BernoulliPhylogenyModel, data::AbstractHLAData, 
                    replacement::Replacement, tree::PhylogeneticTree)

    y = targets(replacement, data)
    m = hla_matrix(data.hla_types)
    ctree = deepcopy(tree)
    annotate!(ctree, data, replacement)

    p = state_probabilities(ctree, TwoStateGTR)
    phylogeny_effect = [p[s]["1"] for s in string.(1:length(y))]
    phylogeny_effect = [min(max(0.01, x), 0.99) for x in phylogeny_effect]

    stan_input = Dict("y" => collect(skipmissing(y)), "hla_matrix" => m[.!ismissing.(y), :],
                      "n_entries" => length(collect(skipmissing(y))), 
                      "n_alleles" => size(m)[2],
                      "phylogeny_effect" => phylogeny_effect[.!ismissing.(y)])
end