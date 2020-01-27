function binned_intervals(
    theta::AbstractVector{<: Real}, y::AbstractVector{<: Bool};
    bins::Int = 10, lower::Real = 0.025, upper = 0.975
)   
    length(theta) == length(y) || 
        throw(DimensionMismatch("theta and y must have same length."))
    
    sort_idx = sortperm(theta)
    theta = theta[sort_idx]
    y = y[sort_idx]

    df = DataFrame(
        expected = Float64[], observed = Float64[], 
        lower = Float64[], upper = Float64[]
    )

    for partition in Iterators.partition(1:length(theta), div(length(theta), bins))
        y_p = y[partition]
        theta_p = theta[partition]
        d = Beta(sum(y_p) + 1, length(y_p) - sum(y_p) + 1)

        interval = invlogcdf.(d, (log(lower), log(upper)))

        push!(df, (mean(theta_p), mean(y_p), interval...))
    end

    return df
end