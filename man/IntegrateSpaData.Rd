\name{IntegrateSpaData}
\alias{IntegrateSpaData}
\title{Integrate multiple SRT data}
\description{
  Integrate multiple SRT data based on the PRECASTObj by PRECAST model fitting.
}
\usage{
  IntegrateSpaData(PRECASTObj, species="Human", 
                 custom_housekeep=NULL, covariates_use=NULL,
                 seuList=NULL, subsample_rate=1, sample_seed=1)
}
\arguments{
  \item{PRECASTObj}{a PRECASTObj object after finishing the PRECAST model fitting and model selection.}
  \item{species}{an optional string, one of 'Human', 'Mouse' and 'Unknown', specify the species of the SRT data to help choose the housekeeping genes. 'Unknown' means only using the PRECAST results reconstruct the alligned gene expression.}
  \item{custom_housekeep}{user-specified housekeeping genes.}
  \item{covariates_use}{a string vector, the colnames  in `PRECASTObj@seulist[[1]]@meta.data`, representing other biological covariates to considered when removing batch effects.  This is achieved by adding additional covariates for biological conditions in the regression, such as case or control. Default as `NULL`, denoting no other covariates to be considered.}
    \item{seuList}{an optional Seurat list object, `seuList` plays a crucial role in the integration process. If `seuList` is set to `NULL` and `PRECASTObj@seuList` is not NULL, then `seuList` will adopt the values of `PRECASTObj@seuList`. Subsequently, the genes within `seuList` will be utilized for integration. Conversely, if `seuList` is not NULL, the integration will directly employ the genes specified within `seuList`. In the event that both `seuList` and `PRECASTObj@seuList` are set to NULL, integration will proceed using the genes outlined in `PRECASTObj@seulist`, i.e., the variable genes. To preserve the `seuList` not NULL in `PRECASTObj@seuList`, user can set `rawData.preserve=TRUE` when running `CreatePRECASTObject`. This parameter empowers users to integrate the entire set of genes in `seuList` when implementing the integration, as opposed to exclusively considering the variable genes within `PRECASTObj@seuList`.}
  \item{subsample_rate}{an optional real number ranging from zero to one, this parameter specifies the subsampling rate during integration to enhance computational efficiency, default as 1 (without subsampling).}
  \item{sample_seed}{an optional integer, with a default value of 1, serves to designate the random seed when `subsample_rate` is set to a value less than one, ensuring reproducibility in the sampling process.}
}
\details{
  Nothing
}
\value{
  Return a Seurat object by integrating all SRT data batches into a SRT data, where the column "batch" in the meta.data represents the batch ID, and the column "cluster" represents the clusters obtained by PRECAST.
}

\references{
  \href{https://www.nature.com/articles/s41467-023-35947-w}{Wei Liu, Liao, X., Luo, Z. et al, Jin Liu* (2023). Probabilistic embedding, clustering, and alignment for integrating spatial transcriptomics data with PRECAST. Nature Communications, 14, 296}

  Gagnon-Bartsch, J. A., Jacob, L., & Speed, T. P. (2013). Removing unwanted variation from high dimensional data with negative controls. Berkeley: Tech Reports from Dep Stat Univ California, 1-112.
}


\author{
  Wei Liu
}
\note{
  nothing
}
\seealso{
  None
}
\examples{
  data(PRECASTObj)
  PRECASTObj <- SelectModel(PRECASTObj)
  seuInt <- IntegrateSpaData(PRECASTObj, species='unknown')
}