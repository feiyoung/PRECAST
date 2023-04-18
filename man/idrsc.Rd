\name{ICM.EM}
\alias{ICM.EM}
\title{ICM-EM algorithm implementation}
\description{
  ICM-EM algorithm for fitting PRECAST model
}
\usage{
  ICM.EM(XList, q, K, AdjList=NULL,  Adjlist_car=NULL, posList = NULL, 
      platform = "ST", beta_grid=seq(0.2,4, by=0.2),maxIter_ICM=6,
      maxIter=20, epsLogLik=1e-5, verbose=TRUE,mix_prop_heter=TRUE, 
      Sigma_equal=FALSE, Sigma_diag=TRUE,error_heter=TRUE, Sp2=TRUE,
      wpca_int=FALSE, int.model='EEE', seed=1,coreNum = 1, coreNum_int=coreNum)
}
\arguments{
 \item{XList}{an M-length list consisting of multiple matrices with class \code{dgCMatrix} or \code{matrix} that specify the log-normalization gene expression matrix for each data sample used for PRECAST model.}
  \item{q}{a positive integer, specify the number of latent features to be extracted, default as 15.}
  \item{K}{a positive  integer allowing scalar or vector, specify the number of clusters in model fitting.}
 \item{AdjList}{an M-length list of sparse matrices with class \code{dgCMatrix}, specify the adjacency matrix used for Potts model in PRECAST. We provide this interface for those users  who would like to define the adjacency matrix by their own.}
 \item{Adjlist_car}{an M-length list of sparse matrices with class \code{dgCMatrix}, specify the adjacency matrix used for CAR model in PRECAST, default as AdjList in the Potts model. We provide this interface for those users  who would like to use the different adjacency matrix in CAR model.}
 \item{posList}{an M-length list composed by spatial coordinate matrix for each data sample.}
 \item{platform}{a string, specify the platform of the provided data, default as "Visium".
  There are many platforms to be supported, including ("Visuim", "ST", "SeqFISH", 'merFISH', 'slide-seqv2', 'seqscope', "HDST"). If AdjList is not given, the The platform helps to calculate the adjacency matrix by defining the neighbors.} 
  \item{beta_grid}{an optional vector of positive value, the candidate set of the smoothing parameter to be searched by the grid-search optimization approach.}
  \item{maxIter_ICM}{an optional positive value, represents the maximum iterations of ICM.}
  \item{maxIter}{an optional positive value, represents the maximum iterations of EM.}
  \item{epsLogLik}{an optional positive vlaue, tolerance vlaue of relative variation rate of the observed pseudo log-loglikelihood value, defualt as '1e-5'.}
  \item{verbose}{an optional logical value, whether output the information of the ICM-EM algorithm.}
  \item{mix_prop_heter}{an optional logical value, specify whether betar are distict, default as \code{TRUE}.}
  \item{Sigma_equal}{an optional logical value, specify whether Sigmaks are equal, default as FALSE.}
  \item{Sigma_diag}{an optional logical value, specify whether Sigmaks are diagonal matrices, default as \code{TRUE}.}
  \item{error_heter}{an optional logical value, whether use the heterogenous error for DR-SC model, default as \code{TRUE}. If \code{error_heter=FALSE}, then the homogenuous error is used for probabilistic PCA model in PRECAST.}
  \item{Sp2}{an optional logical value, whether add the ICAR model component in the model, default as TRUE. We provide this interface for those users  who don't want to include the  ICAR model.}
  \item{wpca_int}{an optional logical value, means whether use the weighted PCA to obtain the initial values of loadings and other paramters, default as \code{FALSE} which means the ordinary PCA  is used.}
  \item{int.model}{an optional string, specify which Gaussian mixture model is used in evaluting the initial values for PRECAST, default as "EEE"; and see \code{Mclust} for more models' names.}
  \item{seed}{an optional integer, the random seed in fitting PRECAST model.}
  \item{coreNum}{an optional positive integer, means the number of thread used in parallel computating.}
  \item{coreNum_int}{an optional positive integer, means the number of cores used in parallel computation for initial values when \code{K} is a vector, default as same as \code{coreNum}.}
}
\details{
  Nothing
}
\value{
  ICM.EM returns a \link{list} with class "SeqKiDRSC_Object" with the number of components equal to the length of \code{K}, where each component includes the model fitting results for one number of cluster and is a list consisting of following components:
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
  ## For illustration, maxIter is set to 4
  resList <- ICM.EM(XList,AdjList = AdjList, maxIter=4,
                   q=q, K=K, verbose=TRUE)
}
}