@userplot Calibration_Plot
@userplot Phylogeny_Calibration
@userplot Proportion_Plot

@recipe function f(c::Phylogeny_Calibration)
    legend := false
    grid --> false
    seriescolor --> "#B2001D"
    markerstrokecolor --> "#B2001D"
    xguide --> "observed event percentage"
    yguide --> "bin midpoint"
    formatter := x -> string(Int(round(x * 100))) * "%"

    (c, c.args...)
end

@recipe function f(
    ::Phylogeny_Calibration, 
    result::Union{HLAModelResult{3}, HLAModelResult{4}}
)   
    N = result.sf.data["N"]
    rs = result.sf.data["rs"]
    idx = result.sf.data["idx"]

    theta = [result.sf.data["phy"][rs[i], idx[i]] for i in 1:N]
    y = result.sf.data["y"]

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
    xguide --> "observed event percentage"
    yguide --> "bin midpoint"
    formatter := x -> string(Int(round(x * 100))) * "%"
    
    (c, c.args...)
end

@recipe function f(c::Calibration_Plot, result::HLAModelResultIO)
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
    N = result.sf.data["N"]
    indices = 1:N
    posterior = StanInterface.extract(result.sf)
    theta = Vector{Float64}(undef, length(indices))
    y = Vector{Bool}(undef, length(indices))

    for (i, idx) in enumerate(indices)
        r = result.sf.data["rs"][idx]
        n = result.sf.data["idx"][idx]
        theta[i] = mean(theta_i(result.sf, posterior, i))
        y[i] = result.sf.data["y"][idx]
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

@recipe function f(c::Proportion_Plot)
    legend := false
    grid --> false
    seriescolor --> "#B2001D"
    markerstrokecolor --> "#B2001D"
    linecolor --> "#B2001D"
    xguide --> "proportion of sequences with replacement"
    yguide --> "expected proportion"
    formatter := x -> string(Int(round(x * 100))) * "%"

    (c, c.args...)
end

@recipe function f(::Proportion_Plot, result::HLAModelResult)
    sf = stanfit(result)
    d = sf.data
    posterior = extract(sf)
    R = d["R"]
    y = d["y"]
    rs = d["rs"]

    df = DataFrame(
        expected = Float64[], observed = Float64[], 
        lower = Float64[], upper = Float64[]
    )

    theta = hcat([theta_i(sf, posterior, i) for i in 1:d["N"]]...)

    for r in 1:R
        idx = findall(x -> x == r, rs)
        observed = mean(y[idx])
        expected = map(x -> mean(rand.(Bernoulli.(x[idx]))), eachrow(theta))

        interval = quantile(expected, (0.025, 0.975))

        push!(df, [mean(expected), observed, interval...])
    end

    # diagonal line
    @series begin
        seriestype := :path
        linecolor := "black"
        linestyle := :dash
        [0, 1], [0, 1]
    end

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
    theta, y;
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

    cutpoints = range(0, 1, length = bins)

    for (lo, hi) in zip(cutpoints[1:end-1], cutpoints[2:end])
        idx = findall(x -> lo <= x <= hi, theta)
        y_p = y[idx]

        push!(df, [mean(y_p), lo + (hi - lo) / 2, lo, hi])
    end

    return df
end

function check_calibration_arguments(x, y)
    if length(x) != length(y)
        error("x ($(length(x))) and y ($(length(y))) must be of same length.")
    elseif !all(0 .<= x .<= 1)
        error("elements of x must be between 0 and 1.")
    end
end