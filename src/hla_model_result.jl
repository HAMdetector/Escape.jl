function replacement_summary(result::HLAModelResult)
    sf = stanfit(result)
    stan_input = Escape.stan_input(result)
    data = hla_data(result)
    posterior = extract(sf)

    alleles = hla_alleles(result)
    replacements = Escape.replacements(result)

    df = DataFrame(
        allele = HLAAllele[],
        position = Int[], 
        replacement = Char[],
        posterior_p = Float64[],
        log_odds_lower_95 = Float64[],
        log_odds_upper_95 = Float64[],
        fisher_p = Float64[],
        in_predicted_epitope = Int[]
    )

    R = stan_input["R"]
    D = stan_input["D"]
    Z = stan_input["Z"]
    
    for r in 1:R
        fisher = Escape.run(Escape.FisherTest(), data, replacements[r])

        for d in 1:D
            beta_hla = posterior["beta_hla.$r.$d"]
            lower, upper = quantile(beta_hla, (0.025, 0.975))
            p_value = fisher.p_values[1][alleles[d]]

            push!(df, 
                [ 
                    alleles[d],
                    position(replacements[r]), 
                    replacement(replacements[r]),
                    mean(beta_hla .> 0),
                    lower,
                    upper,
                    p_value,
                    Z[r, d]
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