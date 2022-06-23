# PRECAST
PRECAST: a probabilistic embedding and clustering with alignment for spatial transcriptomics data integration.

# Installation
"PRECAST" depends on the 'Rcpp' and 'RcppArmadillo' package, which requires appropriate setup of computer. For the users that have set up system properly for compiling C++ files, the following installation command will work.
```{Rmd}
# Method 1: install PRECAST from CRAN
install.packages('PRECAST')



# Method 2: Install PRECAST from Github
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("feiyoung/PRECAST")

# If some dependent packages (such as `scater`) on Bioconductor can not be installed nomrally, use following commands, then run abouve command.
if (!require("BiocManager", quietly = TRUE)) ## install BiocManager
    install.packages("BiocManager")
# install the package on Bioconducter
BiocManager::install(c("scater"))
```


For the users that don't have set up system properly, the following setup on different systems can be referred.
## Setup on Windows system
First, download [Rtools](https://cran.r-project.org/bin/windows/Rtools/); second, add the Rtools directory to the environment variable. Users can follow [here](https://helpdeskgeek.com/windows-10/add-windows-path-environment-variable/#:~:text=Go%20ahead%20and%20click%20on%20the%20Environment%20Variables,you%20have%20to%20decide%20which%20one%20to%20edit) to add Windows PATH Environment Variable.


## Setup on MacOS system


First, install `homebrew` in terminal:
```{Linux}
/bin/zsh -c "$(curl -fsSL https://gitee.com/cunkai/HomebrewCN/raw/master/Homebrew.sh)"
```
Second install Xcode. Installation about Xcode can be referred [here](https://stackoverflow.com/questions/8291146/xcode-installation-on-mac#:~:text=You%20get%20it%20from%20the%20Mac%20App%20Store.,find%20the%20app%2C%20and%20click%20the%20install%20button).

Third, install `Xcode CLT`:
```{Linux}
xcode-select --install
```

Forth, install "gcc" for compiling C++ and Fortran.
```{Linux}
brew install gcc
```

## Setup on Linux  system
If you use conda environment on Linux system and some dependent packages (such as `scater`) can not normally installed, you can search R package at [here](https://anaconda.org/). We take the `scater` package as example, and its search result is [here](https://anaconda.org/bioconda/bioconductor-scater). Then you can install it in conda environment by following command.
```{Linux}

conda install -c bioconda bioconductor-scater
```
For the user not using conda environment, if  dependent packages (such as `scater`) not normally installed are in Bioconductor, then use the following command to install the dependent packages.
```{Linux}
# install BiocManager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
# install the package on Bioconducter
BiocManager::install(c("scater"))
```
If  dependent packages (such as `DR.SC`) not normally installed are in CRAN, then use the following command to install the dependent packages.
```{Linux}
# install the package on CRAN
install.packages("DR.SC")
```

## Other notes

For running big data, users can use the following system command to set the C_stack unlimited in case of `R Error: C stack usage is too close to the limit`.
```{Linux}
ulimit -s unlimited
```

# Demonstration

For an example of typical PRECAST usage, please see our [Package Website](https://feiyoung.github.io/PRECAST/index.html) for a demonstration and overview of the functions included in PRECAST.