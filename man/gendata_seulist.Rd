\name{gendata_seulist}
\alias{gendata_seulist}
\title{Generate simulated data}
\description{
  Generate simulated spatial transcriptomics data with two data batches.
}
\usage{
  gendata_seulist(height1=30, width1=30,height2=height1, width2=width1, 
                           p =100, q=10, K=7,  G=4, beta=1.2, sigma2=1, 
                           alpha=8, seed=1, view=FALSE)
}
\arguments{
 \item{height1, width1, height2, width2}{Height and width of lattice grids for generating spatial coordinates, n=height * width spots for expression matrix data}
 \item{p}{number of genes to generate.}
 \item{q}{number of true latent features to generate gene expression}
 \item{K}{number of clusters (spatial domains/cell types).}
 \item{G}{the number of neighbors. The latter must be one of G = 4 or G = 8, which respectively correspond to a first order and a second order dependency structure. By default, G = 4.}
 \item{sigma2}{Variance of error term in probabilitic PCA model.}
 \item{alpha}{a positive factor of mixture mean values.}
 \item{seed}{random seed for generate data}
 \item{view}{Logical value indicating whether the draw should be printed. Do not display the optional borders.}
}
\details{
  Nothing
}
\value{
 return a list consisting of Seurat objects, where each object is a SRT data batch, and  the metadata of this Seurat object will include two columns with names "row" and "col" which are the spatial coordinates.
}
\references{
  None
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

  ## we generate two spatial transcriptomics data 
  seulist <- gendata_seulist(height1=20, width1=20,p=200, K=4)
  seulist
}
