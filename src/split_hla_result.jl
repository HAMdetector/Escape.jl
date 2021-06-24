struct SplitHLAModelResult
    f::Function
    max_idx::Int

    function SplitHLAModelResult(f, max_idx)
        f(1) isa HLAModelResult || error("f(1) does not return a HLAModelResult.")

        new(f, max_idx)
    end
end

function replacement_summary(
    result::SplitHLAModelResult; 
    fdr::Bool = false,fisher_p::Bool = false
)
    df = DataFrame()
    
    @showprogress for split_idx in 1:result.max_idx
        df = vcat(df, replacement_summary(Base.invokelatest(result.f, split_idx), 
            fisher_p = fisher_p))
    end

    sort!(df, :posterior_p, rev = true)

    return df
end

function Escape.loo(result::SplitHLAModelResult)
    pw = Vector{Loo.PointwiseLoo}()
    draws = 0
    
    @showprogress for split_idx in 1:result.max_idx
        split_result = Base.invokelatest(result.f, split_idx)
        d = stan_input(split_result)
        split_pw = Vector{Loo.PointwiseLoo}(undef, d["N"])
        sf = stanfit(split_result)
        p = extract(sf)
        draws = sf.iter * sf.chains

        @threads for i in 1:d["N"]
            ll = pointwise_loglikelihoods(sf, p, i)
            split_pw[i] = Loo.pointwise_loo(ll)
        end

        append!(pw, split_pw)
    end

    return Loo.LooResult(pw, (draws, length(pw)))
end

@recipe function f(
    ::Calibration_Plot, result::SplitHLAModelResult
)
    y = Vector{Int}()
    theta_pred = Vector{Float64}()

    for split_idx in 1:result.max_idx
        split_result = Base.invokelatest(result.f, split_idx)
        sf = stanfit(split_result)

        N = sf.data["N"]
        rs = sf.data["rs"]
        idx = sf.data["idx"]
        d = sf.data
        posterior = extract(sf)

        split_y = sf.data["y"]
        split_theta_pred = Vector{Float64}(undef, N)

        @threads for i in 1:N
            theta = mean(theta_i(sf, posterior, i))
            split_theta_pred[i] = theta 
        end

        append!(y, split_y)
        append!(theta_pred, split_theta_pred)
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

@recipe function f(
    ::Phylogeny_Calibration, result::SplitHLAModelResult
)

    y = Vector{Int}()
    theta  = Vector{Float64}()

    for split_idx in 1:result.max_idx
        split_result = Base.invokelatest(result.f, split_idx)
        sf = split_result.sf
        N = sf.data["N"]
        rs = sf.data["rs"]
        idx = sf.data["idx"]

        split_theta = [sf.data["phy"][rs[i], idx[i]] for i in 1:N]
        split_y = sf.data["y"]

        append!(y, split_y)
        append!(theta, split_theta)
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