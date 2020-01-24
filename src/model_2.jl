struct Model2 <: HLAModel end

struct Model2Result <: HLAModelResult 
    sf::Stanfit
    alleles::Vector{HLAAllele}
    replacements::Vector{Replacement}
end

function run(model::Model2, data::AbstractHLAData; iter::Int = 1000, 
             chains = 4, warmup::Int = 1000, wp::WorkerPool = WorkerPool(workers()), 
             depth::Int = 1, mincount::Int = 10)

    input = stan_input(model, data, depth = depth, mincount = mincount)
    input["R"] > 0 || error("No replacements found.")

    sf = stan(joinpath(@__DIR__, "..", "data", "stan", "model_2"),
              input, stan_args = "adapt delta=0.97", iter = iter, chains = chains, 
              warmup = warmup, wp = wp, refresh = 1)
    alleles = sort(unique_alleles(data.hla_types, depth = depth))
    r = replacements(data, mincount = mincount)

    return Model2Result(sf, alleles, r)
end

function stan_input(
    model::Model2, data::AbstractHLAData; 
    depth::Int = 1, mincount::Int = 10
)
    X = hla_matrix(data.hla_types; depth = depth)
    r = replacements(data, mincount = mincount)
    N = size(X)[1]
    D = size(X)[2]
    R = length(r)

    ys = Matrix{Int64}(undef, R, N + 1)
    xs = Matrix{Float64}(undef, R, N * D)
    y_mean = Vector{Float64}(undef, R)
    hla_mean = [mean(X[:, i]) for i in 1:D]
    fill!(ys, -1)
    fill!(xs, -1)

    for i in 1:R
        t = targets(r[i], data)
        y_mean[i] = mean(skipmissing(t))
        count = 0
        for j in 1:N
            ys[i, ]
            if !ismissing(t[j])
                count += 1
                ys[i, count + 1] = t[j]
            end
        end

        ys[i, 1] = count
        X_ = X[.!ismissing.(t), :]
        xs[i, 1:D*count] = reshape(X_', :, 1)
    end


    d = Dict("N" => N, "D" => D, "R" => R, "ys" => ys, "xs" => xs,
             "y_mean" => y_mean, "hla_mean" => hla_mean, "p0" => 5)
    
    return d
end

function pointwise_loglikelihoods(result::Model2Result, r::Int, i::Int)
    sf = result.sf
    D = sf.data["D"]
    
    log_lik_i = Vector{Vector{Float64}}(undef, sf.chains)

    for c in 1:sf.chains
        log_lik_i[c] = let
            ll = Float64[]

            for iter in 1:sf.iter
                intercept = sf.result[c]["intercepts.$r"][iter]
                beta_hla = [sf.result[c]["beta_hla.$r.$d"][iter] for d in 1:D]
                x_i = sf.data["xs"][r, :][(1 + (i - 1) * D):(i * D)]

                theta = logistic(intercept + dot(x_i, beta_hla))

                y = sf.data["ys"][r, i + 1]
                push!(ll, logpdf(Bernoulli(theta), y))
            end

            ll
        end
    end

    return log_lik_i
end

function indices(result::Model2Result)
    sf = result.sf
    R = sf.data["R"]
    N = sf.data["ys"][:, 1]

    x = Vector{Tuple{Int, Int}}()
    sizehint!(x, R * maximum(N))

    for r in 1:R
        for n in 1:N[r]
            push!(x, (r, n))
        end
    end

    return x
end