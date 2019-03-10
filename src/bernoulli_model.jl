export BernoulliModel, BernoulliResult

struct BernoulliModel <: HLAModel
    prior::Symbol
    iter::Int
    chains::Int
end

struct BernoulliResult <: HLAModelResult
    sf::Stanfit
    alleles::Vector{HLAAllele}
    replacement::Replacement
end

BernoulliModel(; iter = 2000, chains = 4) = BernoulliModel(:finnish_horseshoe, iter, chains)

function run(model::BernoulliModel, data::AbstractHLAData, replacement::Replacement;
             wp::WorkerPool = WorkerPool(), depth::Int = 1)
    if model.prior == :finnish_horseshoe
        path = joinpath(dirname(@__DIR__), "data", "stan", "bernoulli_hs")
    elseif model.prior == :broad_t
        path = joinpath(dirname(@__DIR__), "data", "stan", "bernoulli")
    elseif model.prior == :narrow_t
        path = joinpath(dirname(@__DIR__), "data", "stan", "bernoulli_narrow_t")
    end

    input = stan_input(model, data, replacement)
    sf = stan(path, input, chains = model.chains, iter = model.iter, wp = wp,
              stan_args = "adapt delta=0.97")
    alleles = sort(unique_alleles(filter(x -> missing ∉ x, data.hla_types), depth = depth))

    return BernoulliResult(sf, alleles, replacement)
end

function stan_input(model::BernoulliModel, data::AbstractHLAData, replacement::Replacement;
                    depth::Int = 1)
    y = targets(replacement, data)
    m = hla_matrix(data.hla_types, depth = depth)

    stan_input = Dict("y" => collect(skipmissing(y)), "hla_matrix" => m[.!ismissing.(y), :],
                      "n_entries" => length(collect(skipmissing(y))), 
                      "n_alleles" => size(m)[2])
end