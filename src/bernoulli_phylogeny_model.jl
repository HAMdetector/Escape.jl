export BernoulliPhylogenyModel, BernoulliPhylogenyResult

struct BernoulliPhylogenyModel <: HLAPhylogenyModel
    prior::Symbol
    iter::Int
    chains::Int

    function BernoulliPhylogenyModel(prior, iter, chains)
        if prior ∉ (:broad_t, :finnish_horseshoe)
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
             replacement::Replacement)
    
    tree = phylogenetic_tree(data)
    run(model, data, replacement, tree)
end

function run(model::BernoulliPhylogenyModel, data::AbstractHLAData, 
             replacement::Replacement, tree::PhylogeneticTree; 
             wp::WorkerPool = WorkerPool())
    if model.prior == :broad_t
        path = joinpath(@__DIR__, "..", "data", "stan", "bernoulli_phylogeny")
    elseif model.prior == :finnish_horseshoe
        path = joinpath(@__DIR__, "..", "data", "stan", "bernoulli_phylogeny_hs")
    end

    input = stan_input(model, data, replacement, tree)
    sf = stan(path, input, chains = model.chains, iter = model.iter, wp = wp,
              stan_args = "adapt delta=0.95")
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

    stan_input = Dict("y" => collect(skipmissing(y)), "hla_matrix" => m[.!ismissing.(y), :],
                      "n_entries" => length(collect(skipmissing(y))), 
                      "n_alleles" => size(m)[2],
                      "phylogeny_effect" => phylogeny_effect[.!ismissing.(y)])
end