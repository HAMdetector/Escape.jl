module Escape

import BioSequences, StatsFuns
using StanInterface, LightGraphs, MetaGraphs, RCall, Suppressor, Optim, DataFrames,
      JuMP, Ipopt, LinearAlgebra, AutoGrad

include("alleles.jl")
include("replacement.jl")
include("dataset.jl")
include("bernoulli_model.jl")
include("phylogenetic_trees.jl")
include("phylogenetic_background.jl")
include("plots.jl")

end
