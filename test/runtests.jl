using Test, Escape, BioSequences, Suppressor, LightGraphs, MetaGraphs,
      StanInterface, Distributed, HypothesisTests, Serialization, DataFrames, Loo,
      Statistics, Distributions, Plots, FASTX, StaticArrays

include("alleles.jl")
include("abstract_hla_model.jl")
include("hla_data.jl")
include("hla_dataset.jl")
include("dataset_rousseau.jl")
include("replacement.jl")
include("phylogenetic_tree.jl")
include("rate_matrix.jl")
include("phylogenetic_background.jl")
include("epitope_prediction.jl")
# include("loo.jl")
include("epitope_map.jl")
include("fisher_test.jl")
include("escape_result.jl")
include("posterior_predictive_checks.jl")