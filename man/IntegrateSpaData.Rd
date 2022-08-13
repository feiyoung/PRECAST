\name{IntegrateSpaData}
\alias{IntegrateSpaData}
\title{Integrate multiple SRT data}
\description{
  Integrate multiple SRT data based on the PRECASTObj by PRECAST model fitting.
}
\usage{
  IntegrateSpaData(PRECASTObj, species="Human", custom_housekeep=NULL)
}
\arguments{
  \item{PRECASTObj}{a PRECASTObj object after finishing the PRECAST model fitting and model selection.}
  \item{species}{an optional string, one of 'Human', 'Mouse' and 'Unknown', specify the species of the SRT data to help choose the housekeeping genes. 'Unknown' means only using the PRECAST results reconstruct the alligned gene expression.}
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