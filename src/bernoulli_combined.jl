struct BernoulliCombined <: HLAModel end
struct BernoulliCombinedTBB <: HLAModel end

struct BernoulliCombinedResult <: HLAModelResult end

function run(model::BernoulliCombined, data::AbstractHLAData; iter::Int = 1000, chains = 4,
             wp::WorkerPool = WorkerPool(workers()), depth::Int = 1)

    input = stan_input(model, data, depth = depth)
    sf = stan(joinpath(@__DIR__, "..", "data", "stan", "bernoulli_combined"),
              input, stan_args = "adapt delta=0.97", iter = iter, chains = chains, wp = wp,
              refresh = 1)
end

function run(model::BernoulliCombinedTBB, data::AbstractHLAData; iter::Int = 1000, 
             chains = 4, wp::WorkerPool = WorkerPool(workers()), depth::Int = 1)

    input = stan_input(model, data, depth = depth)
    sf = stan(joinpath(@__DIR__, "..", "data", "stan", "bernoulli_combined_tbb"),
              input, stan_args = "adapt delta=0.97", iter = iter, chains = chains, wp = wp,
              refresh = 1)
end

function stan_input(model::BernoulliCombined, data::AbstractHLAData; depth::Int = 1)
    m = hla_matrix(data.hla_types; depth = 1)
    r = replacements(data)
    y_matrix = Matrix{Union{Missing, Int64}}(undef, length(data.hla_types), length(r))

    for i in 1:length(r)
        y_matrix[:, i] = targets(r[i], data)
    end

    y = Int[]
    ii = Int[]
    jj = Int[]

    for i in 1:length(data.hla_types)
        for j in 1:length(r)
            if !ismissing(y_matrix[i, j])
                push!(y, y_matrix[i, j])
                push!(ii, i)
                push!(jj, j)
            end
        end
    end

    d = Dict("y" => y, "ii" => ii, "jj" => jj, "N_y" => length(y), 
             "N" => length(data.hla_types), "M" => size(m)[2], "H" => m, "R" => length(r))
             
    return d
end

function stan_input(model::BernoulliCombinedTBB, data::AbstractHLAData; depth::Int = 1)
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