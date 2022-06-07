\name{Add_embed}
\alias{Add_embed}
\title{Add tSNE embeddings for a Seurat object}
\description{
  Run t-SNE dimensionality reduction on selected features. 
}
\usage{
  Add_embed(embed, seu, embed_name='tSNE' , assay = "RNA")
}
\arguments{
  \item{seuInt}{a Seurat object.}
  \item{n_comp}{an optional positive integer, specify the number of features to be extracted.}
  \item{reduction}{an optional string, means which dimensional reduction (e.g. PRECAST, PCA) to use for the tSNE. Default is PCA}
  \item{assay}{Name of assay that that t-SNE is being run on}
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