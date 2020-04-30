module Escape

using Base.Threads
using LinearAlgebra, Statistics, DelimitedFiles, Distributed, Serialization
using StanInterface, LightGraphs, MetaGraphs, Suppressor, DataFrames,
      HypothesisTests, MultipleTesting, BioSequences, FASTX,
      StatsBase, StatsFuns, Distributions, Loo, RecipesBase, Base.Threads,
      StaticArrays, ProgressMeter

include("alleles.jl")
include("abstract_hla_data.jl")
include("abstract_hla_model.jl")
include("replacement.jl")
include("hla_dataset.jl")
include("dataset_rousseau.jl")
include("phylogenetic_tree.jl")
include("hla_data.jl")
include("phylogenetic_background.jl")
include("epitope_prediction.jl")
include("epitope_map.jl")
include("hla_model.jl")
include("hla_model_result.jl")
include("fisher_test.jl")
include("escape_result.jl")
include("loo.jl")
include("posterior_predictive_checks.jl")

end