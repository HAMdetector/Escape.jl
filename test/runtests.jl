using Test
using Escape

using DataFrames
using Distributions
using FASTX
using HypothesisTests
using LightGraphs
using Loo
using MetaGraphs
using Plots
using Serialization
using StanInterface
using StaticArrays
using Statistics
using Suppressor

include("alleles.jl")
include("abstract_hla_model.jl")
include("hla_data.jl")
include("split_hla_data.jl")
include("hla_dataset.jl")
include("dataset_rousseau.jl")
include("replacement.jl")
include("phylogenetic_tree.jl")
include("epitope_prediction.jl")
include("hla_model_result.jl")
include("loo.jl")
include("epitope_map.jl")
include("hla_model.jl")
include("fisher_test.jl")
include("hla_model_result_io.jl")
include("posterior_predictive_checks.jl")