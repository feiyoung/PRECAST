\name{IntegrateSpaData}
\alias{IntegrateSpaData}
\title{Add UMAP embeddings for a Seurat object}
\description{
  Run a weighted PCA dimensionality reduction
}
\usage{
  IntegrateSpaData(PRECASTObj, species="Human", custom_housekeep=NULL)
}
\arguments{
  \item{PRECASTObj}{a PRECASTObj object after finishing the PRECAST model fitting and model selection.}
  \item{species}{an optional string, specify the species of the SRT data to help choose the housekeeping genes.}
  \item{custom_housekeep}{user-specified housekeeping genes.}
}
\details{
  Nothing
}
\value{
  Return a Seurat object by integrating all SRT data batches into a SRT data, where the column "batch" in the meta.data represents the batch ID, and the column "cluster" represents the clusters obtained by PRECAST.
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