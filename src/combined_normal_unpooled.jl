struct CombinedNormalUnpooled <: HLAModel end

struct CombinedNormalUnpooledResult <: HLAModelResult 
    sf::Stanfit
    alleles::Vector{HLAAllele}
    replacements::Vector{Replacement}
end

function run(model::CombinedNormalUnpooled, data::AbstractHLAData; iter::Int = 1000, 
             chains = 4, warmup::Int = 1000, wp::WorkerPool = WorkerPool(workers()), 
             depth::Int = 1)

    input = stan_input(model, data, depth = depth)
    sf = stan(joinpath(@__DIR__, "..", "data", "stan", "combined_normal_unpooled"),
              input, stan_args = "adapt delta=0.97", iter = iter, chains = chains, 
              warmup = warmup, wp = wp, refresh = 1)
    alleles = sort(unique_alleles(filter(x -> missing ∉ x, data.hla_types), depth = depth))
    r = replacements(data)

    return CombinedNormalPooledResult(sf, alleles, r)
end

function stan_input(model::CombinedNormalUnpooled, data::AbstractHLAData; depth::Int = 1)
    return stan_input(CombinedNormalPooled(), data, depth = depth)
end