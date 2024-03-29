# HAMdetector: A Bayesian regression model that integrates information to detect HLA-associated mutations

![](https://github.com/HAMdetector/Escape.jl/actions/workflows/CI.yml/badge.svg)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://hamdetector.github.io/Escape.jl/dev/)
 [![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

Escape.jl (HAMdetector) is a julia package for identifying HLA-associated substitutions based on annotated sequence data (aligned viral sequences paired to host HLA class I data).

HAMdetector makes use of information from epitope prediction via [MHCflurry 2.0](https://github.com/openvax/mhcflurry) and phylogeny (based on [RAxML-NG](https://github.com/amkozlov/raxml-ng)).
The model is fit using [Stan](https://github.com/stan-dev/cmdstan).

See the [documentation](https://hamdetector.github.io/Escape.jl/dev/):

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://HAMdetector.github.io/Escape.jl/dev/)

# Software requirements

- Linux 64-bit
- git
- Julia 1.6 or newer

# Installation

In the Julia REPL, switch to Pkg-mode by typing the `]` character. 
Then add `HAMdetector_registry` to your registry and install the `Escape` package:

```
pkg> registry add https://github.com/HAMdetector/HAMdetector_registry.git
pkg> add Escape
```

To make sure that HAMdetector was installed correctly, a large suite of unit tests can be run with `] test Escape`.

# Reporting Issues and Contributions

If you experience issues during the installation or usage of HAMdetector, please open an issue in this repository.
Contributions in the form of suggestions, improvements to the documentation, or pull requests are extremely welcome.

# License

The HAMdetector source code is distributed under a permissive MIT license (see file LICENSE).
