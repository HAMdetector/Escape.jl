# function Escape.loo(result::EscapeResult)
#     pw = Loo.PointwiseLoo[]
#     size = Tuple{Int, Int}[]

#     for res in result
#         loo = Escape.loo(res)
#         push!(pw, loo.pointwise_loo...)
#         push!(size, loo.size)
#     end

#     return Loo.LooResult(pw, (size[1][1], sum(x -> x[2], size)))
# end

# function Escape.loo(result::HLAModelResult)
#     indices_ = indices(result)
#     pw = Vector{Loo.PointwiseLoo}(undef, length(indices_))
    
#     for i in 1:length(indices_)
#         ll = pointwise_loglikelihoods(result, indices_[i]...)
#         pw[i] = Loo.pointwise_loo(ll)
#     end

#     draws = result.sf.iter * result.sf.chains

#     return Loo.LooResult(pw, (draws, length(pw)))
# end