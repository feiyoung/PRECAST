\name{selectIntFeatures}
\alias{selectIntFeatures}
\title{Select common genes for multiple data batches}
\description{
   selectIntFeatures prioritizes genes based on the number of times they were selected as SVGs in all data baches, and chose the top genes as the input for the analysis.  We broke ties by examining the ranks of the tied genes in each original dataset and taking those with the highest median rank.}
\usage{
  selectIntFeatures(seulist, spaFeatureList, IntFeatures=2000)
}
\arguments{
  \item{seulist}{a list consisting of Seurat objects, where each object is a SRT data batch.}
  \item{spaFeatureList}{an list consisting of SVGs vectors, where each vector is the top SVGs for each SRT data batch obtained by SPARK or SPARK-X.}
  \item{IntFeatures}{the number of common SVGs genes to be chosen.}
}
\details{
  Nothing
}
\value{
  Return a string vector, the selected gene list for integration in PRECAST.
}

\references{
\href{https://www.nature.com/articles/s41467-023-35947-w}{Wei Liu, Liao, X., Luo, Z. et al, Jin Liu* (2023). Probabilistic embedding, clustering, and alignment for integrating spatial transcriptomics data with PRECAST. Nature Communications, 14, 296}
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
