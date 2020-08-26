function add_shrinkage_factors!(result::HLAModelResult{1})
    si = stan_input(result)
    R = si["R"]
    D = si["D"]

    for res_d in result.sf.result
        for r in 1:R, d in 1:D
            res_d["kappa.$r.$d"] = zeros(length(res_d["lp__"]))
        end
    end
end

function add_shrinkage_factors!(result::Union{HLAModelResult{4}, HLAModelResult{5}})
    si = stan_input(result)
    R = si["R"]
    D = si["D"]

    s_js = map(std, eachcol(si["X"]))

    for res_d in result.sf.result
        for r in 1:R
            y = si["y"][si["rs"] .== r]
            pseudo_var = (1 / mean(y)) * (1 / (1 - mean(y)))
            n = length(y)
            tau_0 = (si["p0"] / (D - si["p0"])) * (sqrt(pseudo_var) / sqrt(n))

            for d in 1:D
                tau = @. res_d["aux1_tau.$r"] * sqrt(res_d["aux2_tau.$r"]) * tau_0
                lambda_j = @. res_d["aux1_lambda.$r.$d"] * sqrt(res_d["aux2_lambda.$r.$d"]) *
                    s_js[d] * exp(res_d["b_epi"] * si["Z"][r, d])
                b_j = @. 1 / (1 + n * inv(pseudo_var) * s_js[d]^2 * res_d["c2.$r.$d"])
                a_j_sq = @. n * inv(pseudo_var) * tau^2 * s_js[d]^2
                kappa_j = @. (1 + b_j * a_j_sq * lambda_j^2) / (1 + a_j_sq * lambda_j^2) 
                
                res_d["kappa.$r.$d"] = kappa_j
            end 
        end
    end

    return result
end

function replacement_summary(result::HLAModelResult; fdr::Bool = false)
    add_shrinkage_factors!(result)
    sf = stanfit(result)
    si = stan_input(result)
    data = hla_data(result)
    posterior = extract(sf)

    alleles = hla_alleles(result)
    rs = replacements(result)

    df = DataFrame(
        allele = HLAAllele[],
        position = Int[], 
        replacement = Char[],
        posterior_p = Float64[],
        inclusion_p = Float64[],
        log_odds_lower_95 = Float64[],
        log_odds_upper_95 = Float64[],
        fisher_p = Float64[],
        in_predicted_epitope = Int[],
        n_replacement_total = Int[],
        n_replacement_with_hla = Int[]
    )

    R = si["R"]
    D = si["D"]
    Z = si["Z"]
    fisher = Escape.run(Escape.FisherTest(), data; fdr = fdr)

    for r in 1:R
        fisher_idx = findfirst(x -> x == rs[r], fisher.replacements)

        for d in 1:D
            beta_hla = posterior["beta_hla.$r.$d"]
            lower, upper = quantile(beta_hla, (0.025, 0.975))
            p_value = fisher.p_values[fisher_idx][alleles[d]]
            inclusion_p = median(1 .- posterior["kappa.$r.$d"])
                
            push!(df, 
                [ 
                    alleles[d],
                    position(rs[r]), 
                    replacement(rs[r]),
                    mean(beta_hla .> 0),
                    inclusion_p,
                    lower,
                    upper,
                    p_value,
                    Z[r, d],
                    sum(fisher.counts[fisher_idx][alleles[d]][:, 1]),
                    fisher.counts[fisher_idx][alleles[d]][1, 1]
                ]
            )
        end
    end

    sort!(df, :posterior_p, rev = true)

    return df
end

hla_model(result::HLAModelResult) = result.model
hla_data(result::HLAModelResult) = result.data
stanfit(result::HLAModelResult) = result.sf
stan_input(result::HLAModelResult) = stanfit(result).data
hla_alleles(result::HLAModelResult) = result.alleles
replacements(result::HLAModelResult) = result.replacements