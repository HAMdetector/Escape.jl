using Test, Escape, BioSequences, Suppressor, LightGraphs, MetaGraphs, JLD2, FileIO,
      StanInterface, Distributed, HypothesisTests, Serialization, DataFrames

import LightGraphs.SimpleGraphs
import Escape.EpitopeMap

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
include("epitope_map.jl")
include("bernoulli_model.jl")
include("bernoulli_phylogeny_model.jl")
include("fisher_test.jl")
include("hla_analysis.jl")
include("plots.jl")