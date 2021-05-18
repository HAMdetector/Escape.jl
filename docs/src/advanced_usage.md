# Advanced usage

## Saving intermediary results

When running HAMdetector, two distinct steps are executed behind the scenes:
In a first step, the input data for CmdStan has to be computed. This involves
computing replacement probabilities based on the phylogenetic tree and 
running MHCflurry for epitope prediction.
In a second step, CmdStan is run with the previously computed input data.

To just execute the first step, use the function `Escape.stan_input`, which
accepts the same `alignment`, `tree` and `hla_annotations` arguments as the
`Escape.run_model` function.

This has the advantage that later re-running the Stan part of HAMdetector does not 
require re-computation of the phylogeny and epitope prediction.

An example of calculating the Stan inputs first, saving them to disk and then running the
Stan part of HAMdetector is shown below:

```julia
using Escape, JLD2

input = Escape.stan_input(
    alignment = "/home/user/Desktop/alignment.fasta",
    tree = "/home/user/Desktop/phylogeny.tree",
    hla_annotations = "/home/user/Desktop/hla_annotations.csv"
)

JLD2.save("/home/user/Desktop/stan_input.jld2", Dict("stan_input" => input))

# providing the stan_input argument skips re-runnning phylogeny / epitope prediction:
result = Escape.run_model(
    alignment = "/home/user/Desktop/alignment.fasta",
    tree = "/home/user/Desktop/phylogeny.tree",
    hla_annotations = "/home/user/Desktop/hla_annotations.csv",
    stan_input = input
)
```

## HLAData structures

Keeping track of 3 different inputs file for every analysis can be bothersome.
Therefore, alignment, phylogeny and HLA annotation data can be combined into a 
single structure that contains all required information to run an HAMdetector model.
The underlying data is copied, so a saved HLAData object works even if the files it was
generated with are moved or deleted.

```julia
using Escape

hla_data = Escape.HLAData(
    alignment = "/home/user/Desktop/alignment.fasta",
    tree = "/home/user/Desktop/phylogeny.tree",
    hla_annotations = "/home/user/Desktop/hla_annotations.csv"
)
```

HLAData structs can be supplied as a single positional argument in place of
the alignment, tree and hla_annotations keyword arguments:

```julia
result = Escape.run_model(hla_data)
```

## Running sub-models

By default, the `run_model` function uses the full HAMdetector model that includes
the horseshoe prior, phylogeny information and epitope prediction.
Reduced versions of the full model can be run by supplying an additional model keyword
argument:

```julia
using Escape

result = Escape.run_model(
    alignment = "/home/user/Desktop/alignment.fasta",
    tree = "/home/user/Desktop/phylogeny.tree",
    hla_annotations = "/home/user/Desktop/hla_annotations.csv",
    model = Escape.HLAModel{1}()
)
```

The different versions are:

- model 1: A simple logistic regression model
- model 2: Like model 1, but additionally with horseshoe prior
- model 3: Like model 2, but additionally with phylogeny information included
- model 4: Full HAMdetector model, like model 3 but with epitope prediction included.

Note that the full HAMdetector model included all the previous models as a special case.

## Modifying sampler parameters

CmdStan options can be supplied by additional keyword arguments. By default,
4 chains with 1000 iterations each and a warmup of 300 samples are drawn.
This should be plenty for all use-cases, but can be changed as shown below:

```julia
using Escape

result = Escape.run_model(
    alignment = "/home/user/Desktop/alignment.fasta",
    tree = "/home/user/Desktop/phylogeny.tree",
    hla_annotations = "/home/user/Desktop/hla_annotations.csv",
    chains = 4, iter = 1000, warmup = 300, stan_args = "adapt_delta=0.95"
)
```

## Speeding up computations

Run-time of the phylogeny calculations can be reduced drastically by
either setting the `JULIA_NUM_THREADS` environment variable to the number of available
CPU threads or equivalently by starting julia with `julia -t auto`.

Similarly, you should also make sure to enable multi-threading in CmdStan by specifying
`STAN_THREADS=true` in the make/local file inside the CmdStan directory, as explained
in the [CmdStan manual](https://mc-stan.org/docs/2_26/cmdstan-guide/parallelization.html).