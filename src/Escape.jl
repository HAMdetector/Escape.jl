module Escape

import BioSequences, StatsFuns, StatsBase
using LinearAlgebra, Calculus, Statistics, DelimitedFiles, Distributed
using StanInterface, LightGraphs, MetaGraphs, RCall, Suppressor,
      JuMP, Ipopt, JLD2, FileIO

include("alleles.jl")
include("abstract_hlamodel.jl")
include("dataset.jl")
include("dataset_rousseau.jl")
include("replacement.jl")
include("phylogenetic_tree.jl")
include("rate_matrix.jl")
include("substitution_models.jl")
include("phylogenetic_background.jl")
include("bernoulli_model.jl")
include("bernoulli_phylogeny_model.jl")
include("hla_analysis.jl")
include("plots.jl")

end