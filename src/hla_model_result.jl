function replacement_summary(result::HLAModelResult; fdr::Bool = false)
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
            inclusion_p = mean(beta_hla .> 0)
                
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