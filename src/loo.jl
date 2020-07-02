function Escape.loo(result::HLAModelResultIO)
    pw = Loo.PointwiseLoo[]
    size = Tuple{Int, Int}[]

    for res in result
        loo = Escape.loo(res)
        push!(pw, loo.pointwise_loo...)
        push!(size, loo.size)
    end

    return Loo.LooResult(pw, (size[1][1], sum(x -> x[2], size)))
end

function Escape.loo(result::HLAModelResult)
    d = stan_input(result)
    N = d["N"]
    pw = Vector{Loo.PointwiseLoo}(undef, N)
    
    @showprogress for i in 1:N
        ll = pointwise_loglikelihoods(result.sf, i)
        pw[i] = Loo.pointwise_loo(ll)
    end

    draws = result.sf.iter * result.sf.chains

    return Loo.LooResult(pw, (draws, length(pw)))
end

function pointwise_loglikelihoods(sf::StanInterface.Stanfit, i::Int)
    y = sf.data["y"][i]
    iter = sf.iter
    chains = sf.chains

    ll = Vector{Vector{Float64}}(undef, chains)
    for idx in eachindex(sf.result)
        theta = theta_i(sf, i)
        ll[idx] = y == 1 ? log.(theta) : log.(1 .- theta)
    end

    return ll
end

function theta_i(sf::StanInterface.Stanfit, i::Int)
    p = extract(sf)
    stan_input = sf.data

    rs = stan_input["rs"]
    R = stan_input["R"]
    D = stan_input["D"]
    X = stan_input["X"]
    idx = stan_input["idx"]
    phy = stan_input["phy"]

    r = rs[i]

    tau = @. p["aux1_tau.$r"] * sqrt(p["aux2_tau.$r"]) 
    
    nu = zeros(length(p["lp__"]))

    for d in 1:D
        lambda = @. p["aux1_lambda.$r.$d"] * sqrt(p["aux2_lambda.$r.$d"])
        lambda_tilde = @. sqrt(p["c2.$r.$d"] * lambda^2 / (p["c2.$r.$d"] + tau^2 * lambda^2))
        beta_hla = @. lambda_tilde * tau * p["z_std.$r.$d"]

        nu += X[idx[i], d] .* beta_hla
    end

    theta = @. logistic(p["b0_hla.$r"] + (p["b_phy"] * logit(phy[r, idx[i]])) + nu)

    return theta
end