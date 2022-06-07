# PRECAST
PRECAST: a probabilistic embedding and clustering with alignment for spatial transcriptomics data integration.

# Installation

To install the the packages "PRECAST", firstly, install the 'remotes' package. Besides, "PRECAST" depends on the 'Rcpp' and 'RcppArmadillo' package, which also requires appropriate setting of Rtools and Xcode for Windows and Mac OS/X, respectively.
```{Rmd}
# Install it from github
install.packages("remotes")
remotes::install_github("feiyoung/PRECAST")
```
## Setup on Linux or MacOS system
For running big data, users can use the following system command to set the C_stack unlimited in case of `R Error: C stack usage is too close to the limit`.
```{Linux}
ulimit -s unlimited
```

# Demonstration

For an example of typical PRECAST usage, please see our [Package vignette](https://feiyoung.github.io/PRECAST/index.html) for a demonstration and overview of the functions included in PRECAST.