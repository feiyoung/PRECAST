\name{Add_embed}
\alias{Add_embed}
\title{Add embeddings for a Seurat object}
\description{
  Add embeddings for a Seurat object.
}
\usage{
  Add_embed(embed, seu, embed_name='tSNE' , assay = "RNA")
}
\arguments{
  \item{embed}{an embedding matrix.}
  \item{seu}{a Seurat object.}
  \item{embed_name}{an optional string, the name of embeddings.}
  \item{assay}{Name of assay that that embed is being put}
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
