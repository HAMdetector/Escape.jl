function sample_indices(result::HLAModelResult)
    input = stan_input(result)
    R = input["R"]
    N = input["N"]
    ys = input["ys"]

    indices = Vector{Tuple{Int, Int}}()
    sizehint!(indices, R * N)

    for i in 1:R, j in 1:ys[i, 1]
        push!(indices, (i, j))
    end

    return indices
end

stan_input(result::HLAModelResult) = stanfit(result).data
stanfit(result::HLAModelResult) = result.sf