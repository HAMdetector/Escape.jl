struct NormalPhylogenyPooled <: HLAModel end

struct NormalPhylogenyPooledResult <: HLAModelResult 
    sf::Stanfit
    alleles::Vector{HLAAllele}
    replacements::Vector{Replacement}
    data::Dict
end

function run(model::NormalPhylogenyPooled, ds::HLADataset, dir::String; iter::Int = 1000,
             chains::Int = 4, warmup::Int = 1000, depth::Int = 1)

    for data in ds.data
        input = stan_input(model, data, depth = depth)
        sf = stan(joinpath(@__DIR__, "..", "data", "stan", "normal_phylogeny_pooled"),
                  input, stan_args = "adapt delta=0.97", iter = iter, chains = chains,
                  warmup = warmup, refresh = 1)
        alleles = sort(unique_alleles(filter(x -> missing ∉ x, data.hla_types), 
                       depth = depth))
        r = replacements(data)

        result = NormalPhylogenyPooledResult(sf, alleles, r, input)
        serialize(joinpath(dir, "$(data.name).jls"), result)
    end

    return nothing
end

function run(model::NormalPhylogenyPooled, data::AbstractHLAData; iter::Int = 1000, 
             chains = 4, warmup::Int = 1000, wp::WorkerPool = WorkerPool(workers()), 
             depth::Int = 1)

    input = stan_input(model, data, depth = depth, mincount = 5)
    sf = stan(joinpath(@__DIR__, "..", "data", "stan", "normal_phylogeny_pooled"),
              input, stan_args = "adapt delta=0.97", iter = iter, chains = chains, 
              warmup = warmup, wp = wp, refresh = 1)
    alleles = sort(unique_alleles(filter(x -> missing ∉ x, data.hla_types), depth = depth))
    r = replacements(data)

    return NormalPhylogenyPooledResult(sf, alleles, r, input)
end

function stan_input(
        model::NormalPhylogenyPooled, data::AbstractHLAData; depth::Int = 1, 
        mincount::Int = 5
    )

    m = hla_matrix(data.hla_types; depth = depth)
    r = replacements(data, mincount = mincount)
    N_obs = size(m)[1]
    N_alleles = size(m)[2]
    N_shards = length(r)
    tree = !ismissing(data.tree) ? data.tree : phylogenetic_tree(data)

    y_matrix = Matrix{Union{Missing, Int64}}(undef, N_obs, N_shards)

    for i in 1:N_shards
        y_matrix[:, i] = targets(r[i], data)
    end

    ys = Matrix{Int64}(undef, N_shards, N_obs + 1)
    xs = Matrix{Float64}(undef, N_shards, N_obs * (N_alleles + 1))
    fill!(ys, -10000)
    fill!(xs, -10000)

    for i in 1:N_shards
        t = targets(r[i], data)
        pp = StatsFuns.logit.(Escape.phylogeny_probabilities(r[i], data, tree))

        count = 0
        for j in 1:N_obs
            if !ismissing(t[j])
                count += 1
                ys[i, count + 1] = t[j]
            end
        end
        ys[i, 1] = count

        X = hcat(pp, m[findall(x -> !ismissing(t), t), :])
        xs[i, :] = reshape(X', :, 1)
    end

    d = Dict("N_obs" => N_obs, "N_alleles" => N_alleles, "N_shards" => N_shards,
             "ys" => ys, "xs" => xs)
    
    return d
end