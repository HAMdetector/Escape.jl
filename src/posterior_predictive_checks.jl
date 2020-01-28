@userplot Calibration

# @recipe function f(::HLAModelResult) begin

# end

@recipe function f(c::Calibration)
    x, y = c.args
    check_calibration_arguments(x, y)
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

function check_calibration_arguments(x, y)
    if !(x isa AbstractVector{<: Real})
        error("x must be <: AbstractVector{<: Real}, got $(typeof(x)).")
    elseif !(y isa AbstractVector{<: Bool})
        error("y must be <: AbstractVector{<: Bool}, got $(typeof(y)).")
    elseif length(x) != length(y)
        error("x ($(length(x))) and y ($(length(y))) must be of same length.")
    elseif !all(0 .<= x .<= 1)
        error("elements of x must be between 0 and 1.")
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