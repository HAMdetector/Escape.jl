module Escape

import BioSequences, StatsFuns, StatsBase
using LinearAlgebra, Calculus, Statistics
using StanInterface, LightGraphs, MetaGraphs, RCall, Suppressor,
      JuMP, Ipopt

include("alleles.jl")
include("abstract_hlamodel.jl")
include("dataset.jl")
include("replacement.jl")
include("phylogenetic_tree.jl")
include("rate_matrix.jl")
include("substitution_models.jl")
include("phylogenetic_background.jl")
include("bernoulli_model.jl")
include("bernoulli_phylogeny_model.jl")
include("plots.jl")

end
