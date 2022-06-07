\name{ICM.EM}
\alias{ICM.EM}
\title{ICM-EM algorithm implementation}
\description{
  ICM-EM algorithm for fitting PRECAST model
}
\usage{
  ICM.EM(XList, q, K, AdjList=NULL,  Adjlist_car=NULL, posList = NULL, platform = "ST",
                  beta_grid=seq(0.2,4, by=0.2),
                   maxIter_ICM=6,maxIter=20, epsLogLik=1e-5, verbose=TRUE,
                   mix_prop_heter=TRUE, Sigma_equal=FALSE, Sigma_diag=TRUE,error_heter=TRUE, Sp2=TRUE,
                   wpca_int=FALSE, int.model='EEE', seed=1,coreNum = 1, coreNum_int=coreNum)
}
\arguments{
 \item{XList}{an M-length list consisting of multiple matrices with class \code{dgCMatrix} or \code{matrix} that specify the log-normalization gene expression matrix for each data sample used for iDR-SC model.}
  \item{q}{a positive integer, specify the number of latent features to be extracted, default as 15.}
  \item{K}{a positive  integer allowing scalar or vector, specify the number of clusters in model fitting.}
 \item{AdjList}{an M-length list of sparse matrices with class \code{dgCMatrix}, specify the adjacency matrix used for Potts model in iDR-SC. We provide this interface for those users  who would like to define the adjacency matrix by their own.}
 \item{Adjlist_car}{an M-length list of sparse matrices with class \code{dgCMatrix}, specify the adjacency matrix used for CAR model in iDR-SC, default as AdjList in the Potts model. We provide this interface for those users  who would like to use the different adjacency matrix in CAR model.}
 \item{posList}{an M-length list composed by spatial coordinate matrix for each data sample.}
 \item{platform}{a string, specify the platform of the provided data, default as "Visium".
  There are many platforms to be supported, including ("Visuim", "ST", "SeqFISH", 'merFISH', 'slide-seqv2', 'seqscope', "HDST"). If AdjList is not given, the The platform helps to calculate the adjacency matrix by defining the neighbors.} 
  \item{beta_grid}{an optional vector of positive value, the candidate set of the smoothing parameter to be searched by the grid-search optimization approach.}
  \item{maxIter_ICM}{an optional positive value, represents the maximum iterations of ICM.}
  \item{maxIter}{an optional positive value, represents the maximum iterations of EM.}
  \item{epsLogLik}{an optional positive vlaue, tolerance vlaue of relative variation rate of the observed pseudo log-loglikelihood value, defualt as '1e-5'.}
  \item{verbose}{an optional logical value, whether output the information of the ICM-EM algorithm.}
  \item{mix_prop_heter}{}
  \item{Sigma_equal}{}
  \item{Sigma_diag}{}
  \item{error_heter}{an optional logical value, whether use the heterogenous error for DR-SC model, default as \code{TRUE}. If \code{error_heter=FALSE}, then the homogenuous error is used for probabilistic PCA model in iDR-SC.}
  \item{Sp2}{an optional logical value, whether add the ICAR model component in the model, default as TRUE. We provide this interface for those users  who don't want to include the  ICAR model.}
  \item{wpca.int}{an optional logical value, means whether use the weighted PCA to obtain the initial values of loadings and other paramters, default as \code{FALSE} which means the ordinary PCA  is used.}
  \item{int.model}{}
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
  resList <- ICM.EM(list(dat$X1, dat$X2), posList = NULL, AdjList = AdjList, 
                   q=q, K=K, verbose=TRUE)
}
}