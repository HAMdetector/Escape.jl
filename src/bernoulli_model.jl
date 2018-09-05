export BernoulliModel, BernoulliResult

struct BernoulliModel
    chains::Int
    iter::Int
end

struct BernoulliResult
    sf::Stanfit
    alleles::Vector{HLAAllele}
end

function BernoulliModel()
    return BernoulliModel(4, 2000)
end

function run(model::BernoulliModel, replacement::Replacement, data::HLAData)
    stan_path = joinpath(@__DIR__, "..", "data", "stan", "bernoulli")
    sf = stan(stan_path, stan_input(model, replacement, data), chains = model.chains,
              iter = model.iter)
    alleles = unique_alleles(data.hla_types)

    return BernoulliResult(sf, alleles)
end

function stan_input(model::BernoulliModel, replacement::Replacement, data::HLAData)
    y = targets(replacement, data)
    m = hla_matrix(data.hla_types)

    stan_input = Dict("y" => collect(skipmissing(y)), "hla_matrix" => m[.!ismissing.(y), :],
                      "n_entries" => length(collect(skipmissing(y))), 
                      "n_alleles" => size(m)[2])
end