@userplot Calibration_Plot

function calibration_plot_(result::HLAModelResult)
    indices = Escape.indices(result)
    theta = Vector{Float64}(undef, length(indices))
    y = Vector{Bool}(undef, length(indices))

    for (i, idx) in enumerate(indices)
        theta[i] = mean(Iterators.flatten(thetas(result, idx...)))
        y[i] = result.sf.data["ys"][idx[1], idx[2] + 1]        
    end

    calibration_plot_(theta, y)
end

function calibration_plot_(x::AbstractVector{<: Real}, y::AbstractVector{<: Bool})
    if !(x isa AbstractVector{<: Real})
        error("x must be <: AbstractVector{<: Real}, got $(typeof(x)).")
    elseif !(y isa AbstractVector{<: Bool})
        error("y must be <: AbstractVector{<: Bool}, got $(typeof(y)).")
    elseif length(x) != length(y)
        error("x ($(length(x))) and y ($(length(y))) must be of same length.")
    elseif !all(0 .<= x .<= 1)
        error("elements of x must be between 0 and 1.")
    end

    return x, y
end

@recipe function f(c::Calibration_Plot)
    x, y = calibration_plot_(c.args...)
    df = binned_intervals(x, y)

    legend := false
    grid --> false
    markercolor --> "#B2001D"
    markerstrokecolor --> "#B2001D"
    xlabel --> "observed event percentage"
    ylabel --> "bin midpoint"
    formatter := x -> string(Int(round(x * 100))) * "%"
    
    @series begin
        seriestype := :scatter
        df[!, :expected], df[!, :observed]
    end

    for row in eachrow(df)
        @series begin
            seriestype := :path
            linecolor --> "#B2001D"
            [row[:expected], row[:expected]], [row[:lower], row[:upper]]
        end
    end

    @series begin
        seriestype := :path
        linecolor := "black"
        linestyle := :dash
        [0, 1], [0, 1]
    end
end

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