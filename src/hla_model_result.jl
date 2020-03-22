function sample_indices(result::HLAModelResult)
    ys = stan_input(result)["ys"]

    indices = Vector{Tuple{Int, Int}}()
    sizehint!(indices, sum(ys[:, 1]))

    for i in 1:size(ys)[1]
        N = ys[i, 1]
        for j in 1:N
            push!(indices, (i, j))
        end
    end

    return indices
end

stan_input(result::HLAModelResult) = stanfit(result).data
stanfit(result::HLAModelResult) = result.sf