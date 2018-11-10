export HLAModel, HLAModelResult, HLAPhylogenyModel, classification_accuracy

abstract type HLAModel end
abstract type HLAModelResult end

abstract type HLAPhylogenyModel <: HLAModel end

function classification_accuracy(result::HLAModelResult)
    n_entries = result.sf.data["n_entries"]
    posterior = StanInterface.extract(result.sf)
    
    estimated = [median(posterior["y_rep.$i"]) > 0.5 ? 1 : 0 for i in 1:n_entries]
    observed = result.sf.data["y"]
    
    return Float64(count(estimated .== observed) / n_entries)
end