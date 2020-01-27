module Escape

using Base.Threads
using LinearAlgebra, Calculus, Statistics, DelimitedFiles, Distributed, Serialization
using StanInterface, LightGraphs, MetaGraphs, Suppressor, DataFrames,
      JuMP, Ipopt, HypothesisTests, MultipleTesting, BioSequences,
      StatsBase, StatsFuns, Distributions, Loo

include("alleles.jl")
include("abstract_hla_data.jl")
include("replacement.jl")
include("abstract_hla_model.jl")
include("hla_dataset.jl")
include("loo.jl")
include("dataset_rousseau.jl")
include("phylogenetic_tree.jl")
include("hla_data.jl")
include("rate_matrix.jl")
include("substitution_models.jl")
include("phylogenetic_background.jl")
include("epitope_prediction.jl")
include("epitope_map.jl")
include("model_2.jl")
include("model_3.jl")
include("model_4.jl")
include("fisher_test.jl")
include("escape_result.jl")
include("posterior_predictive_checks.jl")
include("plots.jl")

end