struct Model1 <: HLAModel end

struct Model1Result <: HLAModelResult 
    sf::Stanfit
    alleles::Vector{HLAAllele}
    replacements::Vector{Replacement}
end

function run(
    model::Model1, data::AbstractHLAData;
    mincount::Int = 10,
    depth::Int = 1,
    iter::Int = 1000, 
    chains = 4, 
    warmup::Int = 1000, 
    wp::WorkerPool = WorkerPool(workers())
)
    input = stan_input(model, data, mincount = mincount)
    input["R"] > 0 || error("No replacements found.")

    sf = stan(
        joinpath(@__DIR__, "..", "data", "stan", "model_1"), input,
        stan_args = "adapt delta=0.97",
        iter = iter,
        chains = chains,
        warmup = warmup,
        wp = wp,
        refresh = 1
    )
    alleles = sort(unique_alleles(data.hla_types, depth = depth))
    r = replacements(data, mincount = mincount)[1:20]

    return Model1Result(sf, alleles, r)
end

function stan_input(
    model::Model1, data::AbstractHLAData;
    depth::Int = 1, mincount::Int = 10 
)
    X = hla_matrix(data.hla_types; depth = depth)
    r = replacements(data, mincount = mincount)[1:20]
    N = size(X)[1]
    D = size(X)[2]
    R = length(r)

    ys = Matrix{Int64}(undef, R, N + 1)
    xs = Matrix{Int64}(undef, R, N * D)
    fill!(xs, 0)
    fill!(ys, 0)

    for i in 1:R
        t = targets(r[i], data)
        observations = length(t) - count(ismissing.(t))
        ys[i, 1] = observations
        ys[i, 2:(observations + 1)] = collect(skipmissing(t))

        X_ = X[.!ismissing.(t), :]
        xs[i, 1:(D * observations)] = reshape(X_', :, 1)
    end

    d = Dict("N" => N, "D" => D, "R" => R, "ys" => ys, "xs" => xs)
end

function pointwise_loglikelihoods(result::Model1Result, r::Int, i::Int)
    sf = result.sf
    D = sf.data["D"]
    
    thetas_i = thetas(result, r, i)
    log_lik_i = Vector{Vector{Float64}}(undef, length(thetas_i))

    for c in eachindex(thetas_i)
        log_lik_i[c] = similar(thetas_i[c])

        for n in eachindex(log_lik_i[c])
            y = sf.data["ys"][r, i + 1]
            log_lik_i[c][n] = logpdf(Bernoulli(thetas_i[c][n]), y)
        end
    end

    return log_lik_i
end

function thetas(result::Model1Result, r::Int, i::Int)
    sf = result.sf
    D = sf.data["D"]

    theta_i = Vector{Vector{Float64}}(undef, sf.chains)

    for c in 1:sf.chains
        theta_i[c] = let
            theta = Float64[]

            for iter in 1:sf.iter
                intercept = sf.result[c]["intercepts.$r"][iter]
                beta_hla = [sf.result[c]["beta_hla.$r.$d"][iter] for d in 1:D]
                x_i = sf.data["xs"][r, :][(1 + (i - 1) * D):(i * D)]

                push!(theta, logistic(intercept + dot(x_i, beta_hla)))
            end

            theta
        end
    end

    return theta_i
end

function indices(result::Model1Result)
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