@userplot Calibration_Plot
@userplot Phylogeny_Calibration
@userplot Proportion_Plot

@recipe function f(c::Calibration_Plot)
    legend := false
    grid --> false
    seriescolor --> "#B2001D"
    markerstrokecolor --> "#B2001D"
    xguide --> "observed event percentage"
    yguide --> "expected event percentage"
    formatter := x -> string(Int(round(x * 100))) * "%"
    
    (c, c.args...)
end

@recipe function f(::Calibration_Plot, result::HLAModelResult; bins = 15, samples = 200,
        title = "")
    df = binned_intervals(result, bins = bins, samples = samples)

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
            # linecolor := plotattributes[:linecolor]
            x := [row[:expected], row[:expected]]
            y := [row[:lower], row[:upper]]
            ()
        end
        
    end
end

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

function binned_intervals(result::HLAModelResult; bins = 20, samples = 200)
    stan_input = Escape.stan_input(result)
    
    theta = posterior_predictions(result.sf, stan_input, samples = samples)
    binned_indices = Escape.binned_indices(theta, bins = bins)
    predictions = rand.(Bernoulli.(theta))

    df = DataFrame(
        expected = Float64[], observed = Float64[], 
        lower = Float64[], upper = Float64[]
    )

    for indices in binned_indices
        length(indices) > 0 || continue
        
        expected = map(x -> mean(predictions[x, indices]), 1:samples)
        observed = mean(stan_input["y"][indices])
        lower = minimum(expected)
        upper = maximum(expected)

        push!(df, [mean(expected), observed, lower, upper])
    end

    return df
end

function binned_indices(theta::Matrix; bins = 20)
    cutpoints = range(0, 1, length = bins)
    bins = [Int[] for i in 1:bins]

    for i in 1:size(theta, 2)
        m = median(theta[:, i])

        for (j, (low, high)) in enumerate(zip(cutpoints[1:end-1], cutpoints[2:end]))
            if low <= m <= high
                push!(bins[j], i)
            end
        end
    end

    return bins
end

function posterior_predictions(sf::StanInterface.Stanfit, stan_input::Dict; samples = 200)
    posterior = StanInterface.extract(sf)
    theta = Matrix{Float64}(undef, samples, stan_input["N"])

    p = Progress(size(theta, 2), desc = "theta_i ")

    @threads for i in 1:size(theta, 2)
        theta[:, i] .= Escape.theta_i(sf, posterior, i)[1:samples]
        next!(p)
    end

    return theta
end