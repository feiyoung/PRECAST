\name{IntegrateSpaData}
\alias{IntegrateSpaData}
\title{Integrate multiple SRT data}
\description{
  Integrate multiple SRT data based on the PRECASTObj by PRECAST model fitting.
}
\usage{
  IntegrateSpaData(PRECASTObj, species="Human", 
                 custom_housekeep=NULL, covariates_use=NULL)
}
\arguments{
  \item{PRECASTObj}{a PRECASTObj object after finishing the PRECAST model fitting and model selection.}
  \item{species}{an optional string, one of 'Human', 'Mouse' and 'Unknown', specify the species of the SRT data to help choose the housekeeping genes. 'Unknown' means only using the PRECAST results reconstruct the alligned gene expression.}
  \item{custom_housekeep}{user-specified housekeeping genes.}
  \item{covariates_use}{a string vector, the colnames  in `PRECASTObj@seulist[[1]]@meta.data`, representing other biological covariates to considered when removing batch effects.  This is achieved by adding additional covariates for biological conditions in the regression, such as case or control. Default as `NULL`, denoting no other covariates to be considered.}
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
  PRECASTObj <- selectModel.PRECASTObj(PRECASTObj)
  seuInt <- IntegrateSpaData(PRECASTObj, species='unknown')
}