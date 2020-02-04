@userplot Calibration_Plot
@userplot Phylogeny_Calibration

@recipe function f(c::Phylogeny_Calibration)
    legend := false
    grid --> false
    seriescolor --> "#B2001D"
    markerstrokecolor --> "#B2001D"
    xlabel --> "bin midpoint"
    ylabel --> "observed event percentage"
    formatter := x -> string(Int(round(x * 100))) * "%"

    (c, c.args...)
end

@recipe function f(::Phylogeny_Calibration, result::Union{Model3Result, Model4Result})
    indices = Escape.indices(result)
    D = result.sf.data["D"]
    theta = map(x -> result.sf.data["xs"][x[1], :][1 + (x[2] - 1) * (D + 1)], indices)
    y = map(x -> result.sf.data["ys"][x[1], x[2] + 1], indices)

    @series begin
        seriestype := :calibration
        theta, y
    end

    @series begin
        seriestype := :path
        linecolor := "black"
        linestyle := :dash
        [0, 1], [0, 1]
    end
end

@recipe function f(c::Calibration_Plot)
    legend := false
    grid --> false
    seriescolor --> "#B2001D"
    markerstrokecolor --> "#B2001D"
    xlabel --> "bin midpoint"
    ylabel --> "observed event percentage"
    formatter := x -> string(Int(round(x * 100))) * "%"
    
    (c, c.args...)
end

@recipe function f(c::Calibration_Plot, result::EscapeResult)
    layout --> (ceil(Int, length(result) / 2), 2)

    for (i, res) in enumerate(result)
        @series begin
            title := result.ds.data[i].name
            subplot := i
            (c, res)
        end
    end
    
    plots = ceil(Int, length(result) / 2) * 2
    for i in (length(result) + 1):plots
        @series begin
            subplot := 3
            legend := false
            grid := false
            foreground_color_subplot := :white
            ()
        end
    end
end

@recipe function f(::Calibration_Plot, result::HLAModelResult)
    indices = Escape.indices(result)
    theta = Vector{Float64}(undef, length(indices))
    y = Vector{Bool}(undef, length(indices))

    for (i, idx) in enumerate(indices)
        theta[i] = StatsBase.mode(Iterators.flatten(thetas(result, idx...)))
        y[i] = result.sf.data["ys"][idx[1], idx[2] + 1]        
    end

    @series begin
        seriestype := :calibration
        theta, y
    end

    @series begin
        seriestype := :path
        linecolor := "black"
        linestyle := :dash
        [0, 1], [0, 1]
    end
end

@recipe function f(::Type{Val{:calibration}}, x, y, z)
    bins = plotattributes[:bins] == :auto ? 10 : plotattributes[:bins]
    df = binned_intervals(x, Bool.(y), bins = bins)

    # error bars
    for row in eachrow(df)
        @series begin
            seriestype := :path
            linecolor := plotattributes[:linecolor]
            x := [row[:expected], row[:expected]]
            y := [row[:lower], row[:upper]]
            ()
        end
        
    end

    # scatter points
    seriestype := :scatter
    x := df[!, :expected]
    y := df[!, :observed]
end

function binned_intervals(
    theta::AbstractVector{<:Real}, y::AbstractVector{<:Bool};
    bins::Int = 10, lower::Real = 0.025, upper = 0.975
)   
    check_calibration_arguments(theta, y)
    
    sort_idx = sortperm(theta)
    theta = theta[sort_idx]
    y = y[sort_idx]

    df = DataFrame(
        expected = Float64[], observed = Float64[], 
        lower = Float64[], upper = Float64[]
    )

    binsize = div(length(theta), bins)
    for (i, partition) in enumerate(Iterators.partition(1:length(theta), binsize))
        y_p = y[partition]
        theta_p = theta[partition]
        d = Beta(sum(y_p) + 1, length(y_p) - sum(y_p) + 1)

        interval = invlogcdf.(d, (log(lower), log(upper)))

        if i <= bins
            push!(df, (mean(theta_p), mean(y_p), interval...))
        end
    end

    return df
end

function check_calibration_arguments(x::AbstractVector{<: Real}, y::AbstractVector{<: Bool})
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