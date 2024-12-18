function Escape.loo(result::HLAModelResult)
    d = stan_input(result)
    N = d["N"]
    pw = Vector{Loo.PointwiseLoo}(undef, N)
    p = StanInterface.extract(stanfit(result))

    progress = Progress(N)
    @threads for i in 1:N
        ll = pointwise_loglikelihoods(result.sf, p, i)
        pw[i] = Loo.pointwise_loo(ll)
        next!(progress)
    end
    
    iter = length(p["lp__"])
    chains = length(result.sf.results)

    draws = iter * chains

    return Loo.LooResult(pw, (draws, length(pw)))
end

function pointwise_loglikelihoods(sf::StanInterface.Stanfit, p::Dict, i::Int)
    y = StanInterface.stan_data(sf)["y"][i]

    iter = length(p["lp__"])
    chains = length(sf.results)

    ll = Vector{Vector{Float64}}(undef, chains)
    for idx in eachindex(sf.results)
        theta = theta_i(sf, p, i)
        ll[idx] = y == 1 ? log.(theta) : log.(1 .- theta)
    end

    return ll
end

function theta_i(sf::StanInterface.Stanfit, p::Dict, i::Int)
    stan_input = StanInterface.stan_data(sf)

    rs = stan_input["rs"]
    R = stan_input["R"]
    D = stan_input["D"]
    X = stan_input["X"]
    idx = stan_input["idx"]
    phy = stan_input["phy"]

    r = rs[i]
    
    nu = zeros(length(p["beta_hla.1.1"]))

    for d in 1:D
        beta_hla = p["beta_hla.$r.$d"]

        nu += X[idx[i], d] .* beta_hla
    end

    b_phy = "b_phy.$r" in keys(p) ? p["b_phy.$r"] : zeros(length(p["beta_hla.1.1"]))
    theta = @. logistic(p["b0_hla.$r"] + (b_phy * logit(phy[r, idx[i]])) + nu)

    return theta
end
