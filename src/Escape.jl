module Escape

import BioSequences, StatsFuns, StatsBase
using LinearAlgebra, Calculus, Statistics, DelimitedFiles, Distributed
using StanInterface, LightGraphs, MetaGraphs, Suppressor, DataFrames,
      JuMP, Ipopt, JLD2, FileIO, HypothesisTests, MultipleTesting

include("alleles.jl")
include("abstract_hla_data.jl")
include("replacement.jl")
include("abstract_hla_model.jl")
include("hla_dataset.jl")
include("dataset_rousseau.jl")
include("phylogenetic_tree.jl")
include("hla_data.jl")
include("rate_matrix.jl")
include("substitution_models.jl")
include("phylogenetic_background.jl")
include("bernoulli_model.jl")
include("bernoulli_phylogeny_model.jl")
include("fisher_test.jl")
include("hla_analysis.jl")
include("plots.jl")

end