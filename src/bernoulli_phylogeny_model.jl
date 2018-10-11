export BernoulliPhylogenyModel, BernoulliPhylogenyResult

struct BernoulliPhylogenyModel <: HLAModel
    chains::Int
    iter::Int
end

struct BernoulliPhylogenyResult <: HLAModelResult
    sf::Stanfit
    alleles::Vector{HLAAllele}
end

BernoulliPhylogenyModel(; chains = 4, iter = 2000) = BernoulliPhylogenyModel(chains, iter)

function run(model::BernoulliPhylogenyModel, data::HLAData, replacement::Replacement)
    tree = PhylogeneticTree(data)
    run(model, data, replacement, tree)
end

function run(model::BernoulliPhylogenyModel, data::HLAData, replacement::Replacement, 
             tree::PhylogeneticTree)

    path = joinpath(@__DIR__, "..", "data", "stan", "bernoulli_phylogeny")
    input = stan_input(model, data, replacement, tree)
    sf = stan(path, input, chains = model.chains, iter = model.iter)
    alleles = sort(unique_alleles(data.hla_types))

    return BernoulliPhylogenyResult(sf, alleles)
end


function stan_input(model::BernoulliPhylogenyModel, data::HLAData, replacement::Replacement, 
                    tree::PhylogeneticTree)

    y = targets(data, replacement)
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