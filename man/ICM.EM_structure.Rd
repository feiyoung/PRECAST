\name{ICM.EM_structure}
\alias{ICM.EM_structure}
\title{SRT data integration}
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
  ICM.EM_structure returns a \link{list} with class "SeqKiDRSC_Object" with the number of components equal to the length of \code{K}, where each component includes the model fitting results for one number of cluster and is a list consisting of following components:
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
\href{https://www.biorxiv.org/content/10.1101/2021.12.25.474153v4}{Wei Liu, Xu Liao, Yi Yang, Huazhen Lin, Joe Yeong, Xiang Zhou, Xingjie Shi and Jin Liu. (2022) Joint dimension reduction and clustering analysis for single-cell RNA-seq and spatial transcriptomics data}
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
\dontrun{
  ## we generate the spatial transcriptomics data with lattice neighborhood, i.e. ST platform.
  q <- 10; K <- 4
  dat <- PRECAST:::gendataInte_sp(width1 = 20, height1 = 20,height2=21, width2=21, p=100,q=q, K=K, sigma2 = 9,
                        beta=c(0.5, 1.0))
  AdjList <- lapply(list(dat$pos1, dat$pos2), getAdj_reg, platform='ST')
  message(sapply(AdjList, dim))
  resList <- ICM.EM_structure(list(dat$X1, dat$X2),  AdjList = AdjList, 
                   q=q, K=K)
}
}