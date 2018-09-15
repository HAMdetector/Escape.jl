module Escape

import BioSequences, StatsFuns
using StanInterface, LightGraphs, MetaGraphs, RCall, Suppressor

include("alleles.jl")
include("replacement.jl")
include("dataset.jl")
include("bernoulli_model.jl")
include("phylogenetic_trees.jl")
include("phylogenetic_background.jl")
include("plots.jl")

end
