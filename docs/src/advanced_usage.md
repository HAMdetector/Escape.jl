# Advanced usage

## Running sub-models

By default, the `run_model` function uses the full HAMdetector model that includes
the horseshoe prior, phylogeny information and epitope prediction.
Reduced versions of the full model can be run by supplying an additional positional
argument:

```julia
using Escape

data = HLAData(
    alignment_file = "/home/user/Desktop/alignment.fasta",
    tree_file = "/home/user/Desktop/phylogeny.tree"
)

result = Escape.run_model(HLAModel{4}(), data)
```

The different versions are:

- 1: A simple logistic regression model
- 2: Like model 1, but additionally with horseshoe prior
- 3: Like model 2, but additionally with phylogeny information included
- 4: Full HAMdetector model, like model 3 but with epitope prediction included.

Note that the full HAMdetector model includes all the previous models as a special case.

## Modifying sampler parameters

By default, 4 chains with 1000 iterations each and a warmup of 300 samples are drawn.
This should be plenty for all use-cases, but can be changed as shown below:

```julia
using Escape

data = HLAData(
    alignment_file = "/home/user/Desktop/alignment.fasta",
    tree_file = "/home/user/Desktop/phylogeny.tree"
)

result = Escape.run_model(
    data, chains = 4, iter = 1000, warmup = 300, stan_args = "adapt_delta=0.95"
)
```

For a full list of the CmdStan options, please refer to the [CmdStan sampling parameters](https://mc-stan.org/docs/cmdstan-guide/mcmc-intro.html). It should not be necessary to modify the `stan_args`
argument, it is just listed here for completeness. 

## Using 4-digit HLA alleles

By default, HAMdetector uses HLA alleles with 2-digit accuracy. To use HLA alleles with
4-digit accuracy (if supplied), set the `allele_depth` keyword argument to 2 (default is 1).

```julia
using Escape

data = HLAData(
    alignment_file = "/home/user/Desktop/alignment_with_4_digit_alleles.fasta",
    tree_file = "/home/user/Desktop/phylogeny.tree",
    allele_depth = 2
)

## Speeding up computations

Run-time of the phylogeny calculations can be reduced drastically by
either setting the `JULIA_NUM_THREADS` environment variable to the number of available
CPU threads or equivalently by starting julia with `julia -t auto`.

Similarly, you should also make sure to enable multi-threading in CmdStan by specifying
`STAN_THREADS=true` in the make/local file inside the CmdStan directory, as explained
in the [CmdStan manual](https://mc-stan.org/docs/cmdstan-guide/parallelization.html).