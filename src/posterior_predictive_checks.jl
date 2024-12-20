@userplot Calibration_Plot
@userplot Phylogeny_Calibration

@recipe function f(c::Calibration_Plot)
    legend := false
    grid --> false
    seriescolor --> "#A9A9A9"
    markerstrokecolor --> "#A9A9A9"
    xguide --> "observed event percentage"
    yguide --> "bin midpoint"
    formatter := x -> string(Int(round(x * 100))) * "%"

    (c, c.args...)
end

@recipe function f(
    ::Calibration_Plot, result::HLAModelResult
)
    
    sf = stanfit(result)
    stan_data = StanInterface.stan_data(sf)
    N = stan_data["N"]
    rs = stan_data["rs"]
    idx = stan_data["idx"]
    sf = stanfit(result)

    posterior = StanInterface.extract(sf)
    
    theta_pred = Vector{Float64}(undef, N)
    y = stan_data["y"]

    @threads for i in 1:N
        theta = mean(theta_i(sf, posterior, i))
        theta_pred[i] = theta 
    end

    @series begin
        seriestype := :path
        linecolor := "black"
        linestyle := :dash
        [0, 1], [0, 1]
    end

    @series begin
        seriestype := :calibration
        theta_pred, y
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
    ::Phylogeny_Calibration, result::HLAModelResult
)

if !("phy" ∈ keys(StanInterface.stan_data(result.sf)))
        error("HLAModelResult does not contain phylogeny information.")
    end
    
    stan_data = StanInterface.stan_data(result.sf)
    N = stan_data["N"]
    rs = stan_data["rs"]
    idx = stan_data["idx"]
    
    theta = [stan_data["phy"][rs[i], idx[i]] for i in 1:N]
    y = stan_data["y"]

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

function binned_indices(theta::AbstractMatrix; bins = 20)
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
