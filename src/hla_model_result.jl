function replacement_summary(
    result::HLAModelResult; 
    fdr::Bool = false,
    fisher_p::Bool = false
)
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
        log_odds_lower_95 = Float64[],
        log_odds_upper_95 = Float64[],
        fisher_p = Float64[],
        in_predicted_epitope = Int[],
        n_replacement_total = Int[],
        n_replacement_with_hla = Int[]
    )

    # using length(rs) to circumvent wrong data for blockwise runs
    R = si["R"]
    D = si["D"]
    Z = si["Z"]

    if fisher_p
        fisher = Escape.run(Escape.FisherTest(), data; fdr = fdr)
    end

	for r in 1:length(rs)
        fisher_idx = fisher_p ? findfirst(x -> x == rs[r], fisher.replacements) : 0

        for d in 1:D
            beta_hla = posterior["beta_hla.$r.$d"]
            lower, upper = quantile(beta_hla, (0.025, 0.975))
            p_value = fisher_p ? fisher.p_values[fisher_idx][alleles[d]] : 0
            n_total = fisher_p ? sum(fisher.counts[fisher_idx][alleles[d]][1, :]) : 0
            n_with_hla = fisher_p ? fisher.counts[fisher_idx][alleles[d]][1, 1] : 0
            
            z = Z[r, d]
			# if r > size(Z, 1)
			# 	z = 0
			# else 
			# 	z = Z[r, d]
			# end

            push!(df, 
                [ 
                    alleles[d],
                    position(rs[r]), 
                    replacement(rs[r]),
                    mean(beta_hla .> 0),
                    lower,
                    upper,
                    p_value,
                    z,
                    n_total,
                    n_with_hla
                ]
            )
        end
    end

    sort!(df, :posterior_p, rev = true)

    return df
end

function reduce_size!(result::HLAModelResult)
    reduce_size!(result.sf)

    return result
end

hla_model(result::HLAModelResult) = result.model
hla_data(result::HLAModelResult) = result.data
stanfit(result::HLAModelResult) = result.sf
stan_input(result::HLAModelResult) = stanfit(result).data
hla_alleles(result::HLAModelResult) = result.alleles
replacements(result::HLAModelResult) = result.replacements
