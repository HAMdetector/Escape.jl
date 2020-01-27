using Test, Escape, BioSequences, Suppressor, LightGraphs, MetaGraphs,
      StanInterface, Distributed, HypothesisTests, Serialization, DataFrames, Loo

import LightGraphs.SimpleGraphs
import Escape.EpitopeMap

include("alleles.jl")
include("abstract_hla_model.jl")
include("loo.jl")
include("hla_data.jl")
include("hla_dataset.jl")
include("dataset_rousseau.jl")
include("replacement.jl")
include("phylogenetic_tree.jl")
include("rate_matrix.jl")
include("phylogenetic_background.jl")
include("epitope_prediction.jl")
include("model_2.jl")
include("model_3.jl")
include("model_4.jl")
include("epitope_map.jl")
include("fisher_test.jl")
include("escape_result.jl")
include("plots.jl")