module Escape

using Base.Threads
using LinearAlgebra, Statistics, DelimitedFiles, Distributed, Serialization
using CondaPkg, CSV, StanInterface, Graphs, MetaGraphs, Suppressor, DataFrames,
      HypothesisTests, MultipleTesting, BioSequences, FASTX,
      StatsBase, StatsFuns, StableRNGs, Distributions, Loo, Random, RecipesBase, Base.Threads,
      StaticArrays, ProgressMeter, JLD2, HAMdetector_model_binaries_jll, raxml_ng_jll

include("alleles.jl")
include("abstract_hla_data.jl")
include("replacement.jl")
include("abstract_hla_model.jl")
include("hla_dataset.jl")
include("dataset_rousseau.jl")
include("phylogenetic_tree.jl")
include("hla_data.jl")
include("split_hla_data.jl")
include("phylogenetic_background.jl")
include("epitope_prediction.jl")
include("epitope_map.jl")
include("hla_model.jl")
include("fisher_test.jl")
include("hla_model_result.jl")
include("loo.jl")
include("posterior_predictive_checks.jl")
include("split_hla_result.jl")

function __init__()
    CondaPkg.withenv() do
        run(`mhcflurry-downloads fetch`)
        run(`mhcflurry-downloads fetch models_class1_presentation`)
    end
end

end
