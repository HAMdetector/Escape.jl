function Escape.loo(result::HLAModelResult)
    pw = Loo.PointwiseLoo[]
    sizehint!(pw, length(indices(result)))

    for idx in indices(result)
        ll = pointwise_loglikelihoods(result, idx...)
        push!(pw, Loo.pointwise_loo(ll))
    end

    draws = result.sf.iter * result.sf.chains

    return Loo.LooResult(pw, (draws, length(pw)))
end