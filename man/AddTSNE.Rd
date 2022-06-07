\name{AddTSNE}
\alias{AddTSNE}
\title{Add tSNE embeddings for a Seurat object}
\description{
  Run t-SNE dimensionality reduction on selected features. 
}
\usage{
  AddTSNE(seuInt, n_comp=3, reduction='PRECAST', assay='PRE_CAST', seed=1)
}
\arguments{
  \item{seuInt}{a Seurat object.}
  \item{n_comp}{an optional positive integer, specify the number of features to be extracted.}
  \item{reduction}{an optional string, means which dimensional reduction (e.g. PRECAST, PCA) to use for the tSNE. Default is PCA.}
  \item{assay}{Name of assay that that t-SNE is being run on.}
  \item{seed}{an  optional integer, the random seed to evaluate tSNE.}
}
\details{
  Nothing
}
\value{
  Return a revised Seurat object by adding tSNE reduction object.
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

}