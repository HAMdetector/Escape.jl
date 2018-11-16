using StanInterface

build_binary(joinpath(dirname(@__DIR__), "data", "stan", "bernoulli.stan"),
             joinpath(dirname(@__DIR__), "data", "stan", "bernoulli"))

build_binary(joinpath(dirname(@__DIR__), "data", "stan", "bernoulli_phylogeny.stan"),
             joinpath(dirname(@__DIR__), "data", "stan", "bernoulli_phylogeny"))