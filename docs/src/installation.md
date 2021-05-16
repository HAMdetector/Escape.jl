# Installation

HAMdetector is available on Linux. It depends on [RAxML-NG](https://github.com/amkozlov/raxml-ng)
[MHCflurry](https://github.com/openvax/mhcflurry) and [CmdStan](https://github.com/stan-dev/cmdstan) that need to be installed before running HAMdetector.

## Installing [RAxML-NG](https://github.com/amkozlov/raxml-ng)

For detailed installation instructions please refer to the [RAxML-NG documentation](https://github.com/amkozlov/raxml-ng/wiki/Installation).
In brief, the preferred method of installing RAxML-NG is by 

- downloading the prebuilt binaries from [https://github.com/amkozlov/raxml-ng/releases](https://github.com/amkozlov/raxml-ng/releases). The required version is the one ending in "linux\_x86_64.zip".

- extract the zip archive, move the folder to a path of your liking and make sure that path is added to the system's \$PATH environment variable.

In most cases, this is done by modifying the `~/.profile` file by adding a 
line like `export PATH=$PATH:/directory/containing/the/raxml-ng/binary`, where the path after
the semicolon is modified to point to the directory containing the RAxML-NG binary.

You can test if RAxML-NG is correctly installed by running `raxml-ng --version` from the
command line, which should print the currently installed version number and the list
of developers and contributors.

## Installing [MHCflurry](https://github.com/openvax/mhcflurry)

For detailed installation instructions please refer to the [MHCflurry documentation](https://github.com/openvax/mhcflurry).
In brief, the easiest way to install MHCflurry is by either using pip or conda.

Assuming pip is installed, run

`pip3 install mhcflurry --user` (this command slightly deviates from the official documentation
to avoid common issues when installing one of the dependencies of MHCflurry), followed by

`mhcflurry-downloads fetch`.

You can test if MHCflurry is correctly installed by running `mhcflurry-predict --version`
from the command line, which should print the currently installed version of MHCflurry.

## Installing [CmdStan](https://github.com/stan-dev/cmdstan)

For detailed installation instructions please refer to the [CmdStan documentation](https://mc-stan.org/docs/cmdstan-guide/index.html). 

Assuming a C++ toolchain is already installed,

- download the latest release from [https://github.com/stan-dev/cmdstan/releases(https://github.com/stan-dev/cmdstan/releases). The required version is the one NOT ending in "linux-arm64.tar.gz",
e.g. "cmdstan-2.26.1.tar.gz".

- unpack the tarball to a location of your liking.

- add a file called `local` inside the `make` sub-directory, which contains a single line `STAN_THREADS=true`. Though not strictly required, enabling multi-threading is going
to vastly decrease the run-time of HAMdetector models. 

- from the CmdStan folder, call `make build`. 

- add an environmental variable called `JULIA_CMDSTAN_HOME` pointing to the CmdStan directory.
In most cases, this is done by modifying the `~/.profile` file by adding a line like
`export JULIA_CMDSTAN_HOME=/path/to/the/cmdstan/directory`.

You can test if the environmental variable is correctly set by calling echo $JULIA_CMDSTAN_HOME
from the command line. It should output the location of the CmdStan directory.

## Installing HAMdetector

To install HAMdetector, open a Julia REPL (usually done typing by `julia` from the command line),
switch into package mode by typing `]`, followed by `add https://github.com/HAMdetector/Escape.jl`.

To check if all dependencies are correctly installed and HAMdetector functions as it
is supposed to, run `test Escape`. 