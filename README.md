# Escape.jl

Escape.jl is a package for identifying HLA escape mutations. Based on annotated sequence data (Sanger sequences with HLA types), a Bayesian model is used to detect possible HLA escape mutations.

## Installation

1. install [raxml-ng](https://github.com/amkozlov/raxml-ng). `raxml-ng --version` must work
    from the terminal.
2. install [mhcflurry](https://github.com/openvax/mhcflurry).
    `run(`mhcflurry-predict-scan -h`)` must work from the Julia REPL.
3. clone the repository: `git clone --depth 1 gogs@gogs.zmb.uni-due.de:habermann/Escape.jl.git Escape`
    into the same directory where `StanInterface` and `Loo` are already there.
4. `activate` the project, `instantiate` it (all in Pkg mode).
5. run `build` (in Pkg mode) to compile the Stan models.
6. hopefully, `test` should run without any errors.

## running models

```
ds = HLADataset("Test");
res = Escape.run(Escape.HLAModel{1}(), ds.data[1], warmup = 300, iter = 300, chains = 4)
summary = Escape.replacement_summary(res)
```

Available HLA models right now are `Escape.HLAModel{1}()` and `Escape.HLAModel{4}()`.
Unfortunately, adding custom HLADatasets is kind of difficult right now.
In essence, a HLADataset is a struct with name::String and data::Vector{HLAData} fields.
HLAData structs have the fields

```
name::String
records::Vector{FASTX.FASTA.Record}
hla_types::Vector{HLAType}
tree::Union{Missing, PhylogeneticTree}
stan_input::Union{Missing, Dict{String, Any}}
```

which you basically have to supply yourself right now. The tree and stan_input fields
are optional, adding them avoids recomputing the phylogeny and epitope prediction though.
 
## running models in parallel

By default, Escape runs chains in parallel if more than 1 worker is available.
To activate additional workers, run

```
using Distributed
addprocs(4)
using Escape
@everywhere using Pkg; @everywhere Pkg.activate("Escape"); @everywhere using Escape
```

