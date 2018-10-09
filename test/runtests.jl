using Test, Escape, BioSequences, Suppressor, LightGraphs, MetaGraphs, RCall

include("rate_matrix.jl")

include("alleles.jl")
include("dataset.jl")
include("replacement.jl")
include("bernoulli_model.jl")
include("bernoulli_phylogeny_model.jl")
include("phylogenetic_tree.jl")

include("plots.jl")
