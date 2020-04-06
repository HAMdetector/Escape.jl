stan_input(result::HLAModelResult) = stanfit(result).data
stanfit(result::HLAModelResult) = result.sf