function Escape.loo(result::EscapeResult)
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
        theta = sf.result[idx]["theta.$i"]
        ll[idx] = y == 1 ? log.(theta) : log.(1 .- theta)
    end

    return ll
end