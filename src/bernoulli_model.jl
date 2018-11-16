export BernoulliModel, BernoulliResult

struct BernoulliModel <: HLAModel
    chains::Int
    iter::Int
end

struct BernoulliResult <: HLAModelResult
    sf::Stanfit
    alleles::Vector{HLAAllele}
end

BernoulliModel(; chains = 4, iter = 2000) = BernoulliModel(chains, iter)

function run(model::BernoulliModel, data::HLAData, replacement::Replacement)
    path = joinpath(dirname(@__DIR__), "data", "stan", "bernoulli")
    input = stan_input(model, data, replacement)
    sf = stan(path, input, chains = model.chains, iter = model.iter)
    alleles = sort(unique_alleles(filter(x -> missing ∉ x, data.hla_types)))

    return BernoulliResult(sf, alleles)
end

function stan_input(model::BernoulliModel, data::HLAData, replacement::Replacement)
    y = targets(replacement, data)
    
    complete_cases = findall(i -> !ismissing(y[i]) && missing ∉ data.hla_types[i], 
        1:length(y))
    y = y[complete_cases]
    m = hla_matrix(data.hla_types[complete_cases])

    stan_input = Dict("y" => y, "hla_matrix" => m,
                      "n_entries" => length(y), 
                      "n_alleles" => size(m)[2])
end