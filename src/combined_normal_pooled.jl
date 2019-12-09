struct CombinedNormalPooled <: HLAModel end

struct CombinedNormalPooledResult <: HLAModelResult 
    sf::Stanfit
    alleles::Vector{HLAAllele}
    replacements::Vector{Replacement}
end

function run(model::CombinedNormalPooled, data::AbstractHLAData; iter::Int = 1000, 
             chains = 4, warmup::Int = 1000, wp::WorkerPool = WorkerPool(workers()), 
             depth::Int = 1)

    input = stan_input(model, data, depth = depth)
    sf = stan(joinpath(@__DIR__, "..", "data", "stan", "bernoulli_normal_pooled"),
              input, stan_args = "adapt delta=0.97", iter = iter, chains = chains, 
              warmup = warmup, wp = wp, refresh = 1)
    alleles = sort(unique_alleles(filter(x -> missing ∉ x, data.hla_types), depth = depth))
    r = replacements(data)

    return CombinedNormalPooledResult(sf, alleles, r)
end

function stan_input(model::CombinedNormalPooled, data::AbstractHLAData; depth::Int = 1)
    m = hla_matrix(data.hla_types; depth = 1)
    r = replacements(data)
    N_obs = size(m)[1]
    N_alleles = size(m)[2]
    N_shards = length(r)

    y_matrix = Matrix{Union{Missing, Int64}}(undef, N_obs, N_shards)

    for i in 1:N_shards
        y_matrix[:, i] = targets(r[i], data)
    end

    ys = Matrix{Int64}(undef, N_shards, N_obs + 1)
    xs = Matrix{Float64}(undef, N_shards, N_obs * N_alleles)
    fill!(ys, -10000)
    fill!(xs, -10000)

    for i in 1:N_shards
        t = targets(r[i], data)
        count = 0
        for j in 1:N_obs
            if !ismissing(t[j])
                count += 1
                ys[i, count + 1] = t[j]
            end
        end
        ys[i, 1] = count

        X = m[findall(x -> !ismissing(t), t), :]
        xs[i, :] = reshape(X', :, 1)
    end


    d = Dict("N_obs" => N_obs, "N_alleles" => N_alleles, "N_shards" => N_shards,
             "ys" => ys, "xs" => xs)
    
    return d
end