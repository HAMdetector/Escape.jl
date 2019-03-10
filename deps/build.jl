using StanInterface

build_binary(joinpath(dirname(@__DIR__), "data", "stan", "bernoulli.stan"),
             joinpath(dirname(@__DIR__), "data", "stan", "bernoulli"))

build_binary(joinpath(dirname(@__DIR__), "data", "stan", "bernoulli_phylogeny.stan"),
             joinpath(dirname(@__DIR__), "data", "stan", "bernoulli_phylogeny"))

build_binary(joinpath(dirname(@__DIR__), "data", "stan", "bernoulli_phylogeny_hs.stan"),
             joinpath(dirname(@__DIR__), "data", "stan", "bernoulli_phylogeny_hs"))

build_binary(joinpath(dirname(@__DIR__), "data", "stan", "bernoulli_hs.stan"),
             joinpath(dirname(@__DIR__), "data", "stan", "bernoulli_hs"))

build_binary(joinpath(dirname(@__DIR__), "data", "stan", "bernoulli_narrow_t.stan"),
             joinpath(dirname(@__DIR__), "data", "stan", "bernoulli_narrow_t"))

build_binary(joinpath(dirname(@__DIR__), "data", "stan", "bernoulli_phylogeny_narrow_t.stan"),
             joinpath(dirname(@__DIR__), "data", "stan", "bernoulli_phylogeny_narrow_t"))