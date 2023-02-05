\name{ICM.EM_structure}
\alias{ICM.EM_structure}
\title{ICM-EM algorithm implementation with organized paramters}
\description{
  Efficient data integration as well as  spatial clustering for multiple spatial transcriptomics data
}
\usage{
  ICM.EM_structure(XList,  K, AdjList, q=15,parameterList=NULL)
}
\arguments{
 \item{XList}{an M-length list consisting of multiple matrices with class \code{dgCMatrix} or \code{matrix} that specify the log-normalization gene expression matrix for each data sample used for PRECAST model.}
  \item{K}{a positive  integer allowing scalar or vector, specify the number of clusters in model fitting.}
 \item{AdjList}{an M-length list of sparse matrices with class \code{dgCMatrix}, specify the adjacency matrix used for Potts model and Intrisic CAR model in PRECAST model. We provide this interface for those users  who would like to define the adjacency matrix by their own.}
  \item{q}{a positive integer, specify the number of latent features to be extracted, default as 15.}
  \item{parameterList}{Other arguments in PRECAST model, it can be set by \link{model_set}.}
}
\details{
  Nothing
}
\value{
  ICM.EM_structure returns a \link{list} with class "SeqK_PRECAST_Object" with the number of components equal to the length of \code{K}, where each component includes the model fitting results for one number of cluster and is a list consisting of following components:
\item{cluster}{an M-length list that includes the inferred class labels for each data sample.}
\item{hZ}{an M-length list that includes the batch corrected low-dimensional embeddings for each data sample.}
\item{hV}{an M-length list that includes the estimate the ICAR component for each sample.}
\item{Rf}{an M-length list that includes the posterior probability of domain clusters for each sample.}
\item{beta}{an M-length vector that includes the estimated smoothing parameters for each sample.}
\item{Mu}{mean vectors of mixtures components.}
\item{Sigma}{covariance matrix of mixtures components.}
\item{W}{estimated loading matrix}
\item{Lam}{estimated variance of errors in probabilistic PCA model}
\item{loglik}{pseudo observed log-likelihood.}
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
\examples{
\donttest{
  ## we generate the spatial transcriptomics data with lattice neighborhood, i.e. ST platform.
  library(Matrix)
  q <- 10; K <- 4
  data(PRECASTObj)
  posList <- lapply(PRECASTObj@seulist, function(x) cbind(x$row, x$col))
  AdjList <- lapply(posList, getAdj_reg, platform='ST')
  XList <- lapply(PRECASTObj@seulist, function(x) t(x[['RNA']]@data))
  XList <- lapply(XList, scale, scale=FALSE)
  parList <- model_set(maxIter=4)
  resList <- ICM.EM_structure(XList,  AdjList = AdjList, 
                   q=q, K=K, parameterList=parList)
}
}