struct Model4 <: HLAModel end

struct Model4Result <: HLAModelResult 
    sf::Stanfit
    alleles::Vector{HLAAllele}
    replacements::Vector{Replacement}
end

function run(model::Model4, data::AbstractHLAData; iter::Int = 1000, 
             chains = 4, warmup::Int = 1000, wp::WorkerPool = WorkerPool(workers()), 
             depth::Int = 1, mincount::Int = 10)

    input = stan_input(model, data, depth = depth, mincount = mincount)
    sf = stan(joinpath(@__DIR__, "..", "data", "stan", "model_4"),
              input, stan_args = "adapt delta=0.97", iter = iter, chains = chains, 
              warmup = warmup, wp = wp, refresh = 1)
    alleles = sort(unique_alleles(data.hla_types, depth = depth))
    r = replacements(data, mincount = mincount)

    return Model4Result(sf, alleles, r)
end

function stan_input(
    model::Model4, data::AbstractHLAData; 
    depth::Int = 1, mincount::Int = 10
)   
    tree = ismissing(data.tree) ? phylogenetic_tree(data) : data.tree
    X = hla_matrix(data.hla_types; depth = depth)
    r = replacements(data, mincount = mincount)
    Z = Escape.epitope_feature_matrix(data)[map(x -> x.position, r), :]
    N = size(X)[1]
    D = size(X)[2]
    R = length(r)

    ys = Matrix{Int64}(undef, R, N + 1)
    xs = Matrix{Float64}(undef, R, N * (D + 1))
    y_mean = Vector{Float64}(undef, R)
    hla_mean = [mean(X[:, i]) for i in 1:D]
    fill!(ys, -1)
    fill!(xs, -1)

    phy = @distributed hcat for i in 1:R
        phylogeny_probabilities(r[i], data, tree)
    end

    for i in 1:R
        t = targets(r[i], data)
        y_mean[i] = mean(skipmissing(t))
        count = 0
        for j in 1:N
            if !ismissing(t[j])
                count += 1
                ys[i, count + 1] = t[j]
            end
        end

        ys[i, 1] = count
        X_ = hcat(phy[.!ismissing.(t), i], X[.!ismissing.(t), :])
        xs[i, 1:((D + 1) * count)] = reshape(X_', :, 1)
    end


    d = Dict("N" => N, "D" => D, "R" => R, "ys" => ys, "xs" => xs,
             "y_mean" => y_mean, "hla_mean" => hla_mean, "p0" => 5,
             "Z" => Z)
    
    return d
end

function pointwise_loglikelihoods(result::Model4Result, r::Int, i::Int)
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

function thetas(result::Model4Result, r::Int, i::Int)
    sf = result.sf
    D = sf.data["D"]
    xs = sf.data["xs"]
    x_i = xs[r, (2 + (i - 1) * (D + 1)):(i * (D + 1))]
    phylogeny = xs[r, (1 + (i - 1) * (D + 1))]

    theta_i = [zeros(sf.iter) for i in 1:sf.chains]
    beta_hla = Matrix{Float64}(undef, sf.iter, D)

    for c in 1:sf.chains
        alphas = sf.result[c]["intercepts.$r"] .+ sf.result[c]["beta.1"] .* logit(phylogeny)
        for d in 1:D
            beta_hla[:, d] = sf.result[c]["beta_hla.$r.$d"]
        end

        for iter in 1:sf.iter
            theta_i[c][iter] = logistic(alphas[iter] + dot(beta_hla[iter,:], x_i))
        end
    end

    return theta_i
end

function indices(result::Model4Result)
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