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
    indices = sample_indices(result)
    pw = Vector{Loo.PointwiseLoo}(undef, length(indices))
    
    for i in 1:length(indices)
        ll = pointwise_loglikelihoods(result.sf, indices[i])
        pw[i] = Loo.pointwise_loo(ll)
    end

    draws = result.sf.iter * result.sf.chains

    return Loo.LooResult(pw, (draws, length(pw)))
end

function pointwise_loglikelihoods(sf::StanInterface.Stanfit, idx::Tuple{Int, Int})
    y = sf.data["ys"][idx[1], idx[2] + 1]
    iter = sf.iter
    chains = sf.chains

    ll = Vector{Vector{Float64}}(undef, chains)
    for i in eachindex(sf.result)
        theta = sf.result[i]["theta.$(idx[1]).$(idx[2])"]
        ll[i] = y == 1 ? log.(theta) : log.(1 .- theta)
    end

    return ll
end