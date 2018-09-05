module Escape

import BioSequences
using StanInterface

include("alleles.jl")
include("replacement.jl")
include("dataset.jl")
include("bernoulli_model.jl")

end
