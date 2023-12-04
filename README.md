# PRECAST

=========================================================================
<!-- badges: start -->

[![](https://www.r-pkg.org/badges/version-ago/PRECAST)](https://cran.r-project.org/package=PRECAST)
[![](https://cranlogs.r-pkg.org/badges/PRECAST?color=orange)](https://cran.r-project.org/package=PRECAST)
[![](https://cranlogs.r-pkg.org/badges/grand-total/PRECAST?color=orange)](https://cran.r-project.org/package=PRECAST)
[![DOI](https://zenodo.org/badge/500674213.svg)](https://zenodo.org/badge/latestdoi/500674213)
<!-- badges: end -->

PRECAST: a probabilistic embedding and clustering with alignment for spatial transcriptomics data integration.

PRECAST  is a package for integrating and analyzing multiple spatially resolved transcriptomics (SRT) datasets, developed by the Jin Liu's lab. It unifies spatial factor analysis simultaneously with spatial clustering and embedding alignment, requiring only partially shared cell/domain clusters across datasets.

Check out our [Nature Communications paper](https://www.nature.com/articles/s41467-023-35947-w) and  our [Package Website](https://feiyoung.github.io/PRECAST/index.html) for a more complete description of the methods and analyses. 

PRECAST can be used to compare and contrast experimental datasets in a variety of contexts, for instance:

* Across experimental batches
* Across individuals
* Across different conditions (i.e., case and control)
* Across datasets with only partially shared cell/domain clusters

Once multiple datasets are integrated, the package provides functionality for further data exploration, 
analysis, and visualization. Users can:

* Identify clusters using all data information
* Extract aligned low-dimensional embeddings across datasets
* Recover comparable gene expression matrices among datasets
* Find significant shared (and dataset-specific) gene markers
* Conditional spatially variational genes analysis
* Compare clusters with previously identified domain/cell types
* Visuzlize extracted embeddings using 3-dim tSNE and UMAP
* Visualize clusters and gene expression using tSNE and UMAP

# Installation
"PRECAST" depends on the 'Rcpp' and 'RcppArmadillo' package, which requires appropriate setup of computer. For the users that have set up system properly for compiling C++ files, the following installation command will work.
```{Rmd}
# Method 1: install PRECAST from CRAN
install.packages('PRECAST')


# For the newest version of PRECAST, users can use method 2 for installation.
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



## Usage
For usage examples and guided walkthroughs, check the `vignettes` directory of the repo. 

* [Single SRT data analysis](https://feiyoung.github.io/PRECAST/articles/PRECAST.DLPFC.html)
* [Toy examples for integrating three batches](https://feiyoung.github.io/PRECAST/articles/PRECAST.Simu.html)
* [Integration across experimental batches: DLPFC](https://feiyoung.github.io/PRECAST/articles/PRECAST.DLPFC4.html)
* [Integration across experimental batches: breast cancer](https://feiyoung.github.io/PRECAST/articles/PRECAST.BreastCancer.html)



For the users that don't have set up system properly, the following setup on different systems can be referred.
## Setup on Windows system
First, download [Rtools](https://cran.r-project.org/bin/windows/Rtools/); second, add the Rtools directory to the environment variable.


## Setup on MacOS system
First, install Xcode. Installation about Xcode can be referred [here](https://stackoverflow.com/questions/8291146/xcode-installation-on-mac).


Second, install "gfortran" for compiling C++ and Fortran at [here](https://github.com/fxcoudert/gfortran-for-macOS).


## Setup on Linux  system
For parallel computation on Linux, users must use the following system command to set the C_stack unlimited in case of the error `R Error: C stack usage is too close to the limit`.
```{Linux}
ulimit -s unlimited
```

If you use conda environment on Linux system and some dependent packages (such as `scater`) can not normally installed, you can search R package at anaconda.org website. We take the `scater` package as example, and its search result is https://anaconda.org/bioconda/bioconductor-scater. Then you can install it in conda environment by following command.
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


# Demonstration

For an example of typical PRECAST usage, please see our [Package Website](https://feiyoung.github.io/PRECAST/index.html) for a demonstration and overview of the functions included in PRECAST.

# NEWs

PRECAST version 1.6.2 (2023-08-02)
* Update the code to ensure compatibility with Seurat V5.

PRECAST version 1.6.1 (2023-05-06)
* Fix the bug, reported by `agelber-ucsd`, in `IntegrateSpaData()` when adusting additional covariates.

* Note: the package in CRAN is not updated.


PRECAST version 1.6 (2023-04-18)

* Fix the bug in `PRECAST:::filter_gene()`.

* Revise the function name `selectModel()` to `SelectModel()`, to avoid the mask when loading DR.SC package.

* Update the tutorials.

PRECAST version 1.5 (2023-03-05)

* Fix the [issue](https://github.com/feiyoung/PRECAST/issues/2) reported by anvaly. Specifically, the assay name "RNA" used in functions `IntegrateSpaData()`  is replaced  by  the default assay using `DefaultAssay` function in Seurat. Fix the typo `human = {intersect((genelist),Mouse_HK_genes$Gene)}` with replacement of `human = {intersect(toupper(genelist), PRECAST::Human_HK_genes$Gene)}`。


* Add two functions for visualization: `chooseColors()` and `drawFigs()`.

PRECAST version 1.3 (2022-10-05)

* Fix the [issue](https://github.com/feiyoung/PRECAST_Analysis/issues/1) reported by Boyi Guo. Specifically, the assay name "RNA" used in functions `CreatePRECASTObject()` and  `PRECAST()` is replaced  by  the default assay using `DefaultAssay` function in Seurat.

* Provide more detailed help file for `CreatePRECASTObject()` function. Users can use `?CreatePRECASTObject` in Rstudio to access the help file.
In detail, seuList is a list  with Seurat object as component, and each Seurat object at least includes the raw expression count matrix, and spatial coordinates in metadata for each data batch, where the spatial coordinates information must be saved in the metadata of Seurat, named "row" and "col" for each data batch. See the help file for more details.

* Add the [data](https://github.com/feiyoung/PRECAST/tree/main/vignettes_data) used in [Package Website](https://feiyoung.github.io/PRECAST/index.html).

* Add the wrapper functions for different common-used R objects in the spatial transcriptomics, such as spatialExperiment, etc; see the functions `spe2seurat`, `spe2seuratList`,`seu2seuList` and their help file in our developed [SRTtools R package](https://github.com/feiyoung/SRTtools). This R package is designed to provide auxiliary functions for the serirs of our developed methods, such as [SC-MEB](https://github.com/Shufeyangyi2015310117/SC.MEB), [DR-SC](https://github.com/feiyoung//DR.SC) and [PRECAST](https://github.com/feiyoung/PRECAST).

