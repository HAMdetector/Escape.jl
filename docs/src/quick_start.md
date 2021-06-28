# Quick start

## Running HAMdetector

HAMdetector requires two things: A set of aligned amino acid sequences that are annotated 
with HLA alleles and a phylogenetic tree. The sequence alignment has to be provided
in FASTA format, the phylogenetic tree in Newick format.

The easiest way of running HAMdetector is by providing the input files in a standardized
way, that is, the sequences in the .fasta file are annotated by an underscore separated
list of HLA alleles and the tips of the phylogenetic tree are labeled by their position
in the multiple sequence alignment.

An example of valid input files is shown below:

**alignment\.fasta**
```txt
>HLA-A*32_HLA-A*03_HLA-B*52_HLA-B*14_HLA-C*16_HLA-C*16
MG-RASVMG-RA*V
>HLA-A*68_HLA-A*31_HLA-B*44_HLA-B*08_HLA-C*07_HLA-C*07
MGARASVMG-RT*V
>HLA-A*24_HLA-A*24_HLA-B*15_HLA-B*46_HLA-C*12_HLA-C*12
MSARASVMGSRSSV
>HLA-A*02_HLA-A*03_HLA-B*58_HLA-B*44_HLA-C*14_HLA-C*06
MSARASVMG-RRSV
>HLA-A*03_HLA-A*31_HLA-B*55_HLA-B*40_HLA-C*16_HLA-C*03
MSARASVMG-RMSV
```

**phylogeny.tree**
```txt
((2:0.033310,1:0.118635):33.615353,(3:0.107192,5:0.857526):0.102147,4:0.247706);
```

Assuming these files were saved as `/home/user/Desktop/alignment.fasta` and
`/home/user/Desktop/phylogeny.tree`, HAMdetector can be run with the following commands:

```julia
using Escape

data = HLAData(
    alignment_file = "/home/user/Desktop/alignment.fasta",
    tree_file = "/home/user/Desktop/phylogeny.tree"
)

result = run_model(data)

summary = replacement_summary(result)
```

In a first step, the input data are converted into a suitable format to run the
HAMdetector model. This step includes calculating phylogeny information and running
MHCFlurry for epitope prediction.

In a second step, the HAMdetector model is run on the input data.

The `replacement_summary` function then computes a table of all observed replacements, sorted
by decreasing estimated probability of being HLA allele associated. The column *posterior_p*
denotes the posterior probability of the HLA regression coefficient being larger than 0,
conditioned on the model and the observed data.
A posterior probability of 1 denotes that the model is certain that a given replacement is an HLA-associated mutation, whereas a posterior probability of 0.5 means that there is
no evidence of the replacement being HLA-associated (0.5 because a regression coefficient centered around 0 allocates equal probability to negative and positive values). 

The returned structure is a `DataFrame` object and elements can be accessed using
standard indexing, e.g. `summary[1:10, :posterior_p]` returns the posterior probabilities
of the first ten entries.

DataFrames can be saved as .csv files using the CSV package:

```julia
using CSV

CSV.write("/home/user/Desktop/summary.csv", summary)
```

The `HLAData` and `run_model` functions have an optional keyword argument `save_file`
that can be used to save the converted input data and the result object.

```julia
data = HLAData(
    alignment_file = "/home/user/Desktop/alignment.fasta",
    tree_file = "/home/user/Desktop/phylogeny.tree",
    save_file = "/home/user/Desktop/input_data.jld2"
)

result = run_model(data, save_file = "/home/user/Desktop/model_result.jld2")

summary = replacement_summary(result)
```

The file extension ".jld2" is used because HAMdetector uses the [JLD2](https://github.com/JuliaIO/JLD2.jl) package for saving and loading data.

To restore saved results from disk, use the `load_data` and `load_result` functions:

```julia
using Escape

result = load_result("/home/user/Desktop/model_result.jld2")
data = load_data("/home/user/Desktop/input_data.jld2")
```

## Providing an HLA annotation file

For some datasets, sequences are not annotated with HLA information directly,
but with a unique sequence identifier.
For these cases, the `HLAData` function accepts an additional argument to provide
a file that matches these sequence identifiers with the HLA annotation.
The HLA annotation file has to be a comma-, tab- or semicolon-separated list
of values with the columns "identifier", "A1", "A2", "B1", "B2", "C1", "C2".
An example of a valid annotation file is shown below:

**hla_annotations.csv**
```text
idenfier,A1,A2,B1,B2,C1,C2
seq1,HLA-A*74,HLA-A*24,HLA-B*44,HLA-B*40,HLA-C*04,HLA-C*02
seq2,HLA-A*11,HLA-A*11,HLA-B*51,HLA-B*07,HLA-C*01,HLA-C*01
seq3,HLA-A*11,HLA-A*11,HLA-B*15,HLA-B*35,HLA-C*05,HLA-C*07
seq4,HLA-A*24,HLA-A*23,HLA-B*54,HLA-B*41,HLA-C*06,HLA-C*05
seq5,HLA-A*31,HLA-A*01,HLA-B*57,HLA-B*35,HLA-C*01,HLA-C*06
```

Note that, as with providing HLA annotation directly in the .fasta file, the
order of the individual HLA alleles does not matter, e.g. a sequence annotated with
HLA-B\*18 in the B1 column and HLA-B\*51 in the B2 column is treated the same
as an otherwise identical sequence annotated with HLA-B\*51 in the B1 column and
HLA-B\*18 in the B2 column.

HLA alleles do not have to be provided according to the standard HLA nomenclature,
as HAMdetector is able to parse a broad range of commonly used HLA formats, also see section
[Parsing HLA section](##parsing-hla-alleles).

```julia
using Escape

data = HLAData(
    alignment_file = "/home/user/Desktop/alignment_with_ids.fasta",
    tree_file = "/home/user/Desktop/phylogeny.tree",
    hla_annotation_file = "/home/user/Desktop/hla_annotations.csv"
)
```

If an HLA annotation file is provided, the tips of the phylogenetic tree may either
be labeled by the position of the corresponding sequence in the multiple sequence alignment -
as explained previously- or by the unique sequence identifier.

## Convergence diagnostics

HAMdetector uses CmdStan for model fitting, which in turn uses a variant of
Hamiltonian Monte Carlo sampling to generate samples from the typical set. 
To check that these samples faithfully reflect the posterior distribution,
HAMdetector uses CmdStan's diagnostic tool. The output can be retrieved with:

```julia
using Escape

Escape.diagnostics(result)
```

The output is not nicely formated yet, but a successful model run should end in "Processing complete, no problems detected".

## Posterior predictive checks

Bayesian statistics allows for a form of model checking called [posterior predictive checks](https://arxiv.org/abs/1709.01449).
A simple form of posterior predictive checks are calibration plots, which can be generated with:

```julia
using Escape, Plots

Escape.calibration_plot(result)
```
 
For a well calibrated model, the vertical bars should span the diagonal line.
Strong deviations from this suggest that something went wrong.
A similar plot that specifically checks calibration of the internally calculated
phylogeny can be generated with `Escape.phylogeny_calibration(result)`.

Based on the quality of the provided phylogenetic tree, the calibration plot might
show stronger deviations. This is usually not an issue because HAMdetector is able to handle noisy phylogenetic inputs.

## Parsing HLA alleles

For many datasets, HLA annotations do not follow the [official HLA nomenclature](http://hla.alleles.org/nomenclature/naming.html). HAMdetector allows parsing of a wide range of HLA allele notations and tries to
be as flexible as possible.
To check if HAMdetector is able to parse HLA allele annotations in your dataset, use
the `parse_allele` function:

```julia
using Escape

Escape.parse_allele("B51")
Escape.parse_allele("B5101")
Escape.parse_allele("A*0103")
Escape.parse_allele("Cw*0302")
Escape.parse_allele("HLA-B51")
```

Please note that HAMdetector right now supports either 2-digit or 4-digit HLA alleles,
but no mixture of both. This limitation can in principle be solved by applying a technique called
[partial pooling](http://www.stat.columbia.edu/~gelman/research/published/multi2.pdf), 
which is going to be implemented in future release of HAMdetector.

Missing HLA alleles are supported by providing a blank HLA entry. Internally,
missing HLA alleles are not imputed, but treated as 0 entries in the HLA matrix.
This means that any effect that cannot be explained by one of the present HLA alleles
is accounted for by the error term (the intercept), which is a reasonable default if
HLA alleles are missing at random.