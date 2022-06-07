\name{model_set}
\alias{model_set}
\title{PRECAST model setting}
\description{
  Set the PRECAST model structure and paramters in the algorithm.
}
\usage{
  model_set(Sigma_equal=FALSE, Sigma_diag=TRUE,mix_prop_heter=TRUE,
                      error_heter=TRUE, Sp2=TRUE, wpca_int=FALSE,int.model='EEE',
                      coreNum = 1, coreNum_int=coreNum,
                      beta_grid=seq(0.2,4, by=0.2),
                      maxIter_ICM=6,maxIter=20, epsLogLik=1e-5, verbose=TRUE, seed=1)
}
\arguments{
  \item{Sigma_equal}{an optional logical value, specify whether Sigmaks are equal, default as FALSE.}
  \item{Sigma_diag}{an optional logical value, specify whether Sigmaks are diagonal matrices, default as \code{TRUE}.}
  \item{mix_prop_heter}{an optional logical value, specify whether betar are distict, default as \code{TRUE}.}
  \item{error_heter}{an optional logical value, whether use the heterogenous error i.e. lambdarj != lambdark for each sample r, default as \code{TRUE}. If \code{error_heter=FALSE}, then the homogenuous error is used for probabilistic PCA model.}
  \item{Sp2}{an optional logical value, whether add the ICAR model component in the model, default as TRUE. We provide this interface for those users  who don't want to include the  ICAR model.}
  \item{wpca.int}{an optional logical value, means whether use the weighted PCA to obtain the initial values of loadings and other paramters, default as \code{FALSE} which means the ordinary PCA  is used.}
  \item{int.model}{an optional string, specify which Gaussian mixture model is used in evaluting the initial values for PRECAST, default as "EEE"; and see \code{Mclust} for more models' names.}
  \item{coreNum}{an optional positive integer, means the number of thread used in parallel computating.}
  \item{coreNum_int}{an optional positive integer, means the number of cores used in parallel computation for initial values when \code{K} is a vector, default as same as \code{coreNum}.}
   \item{beta_grid}{an optional vector of positive value, the candidate set of the smoothing parameter to be searched by the grid-search optimization approach.}
  \item{maxIter_ICM}{an optional positive value, represents the maximum iterations of ICM.}
  \item{maxIter}{an optional positive value, represents the maximum iterations of EM.}
  \item{epsLogLik}{an optional positive vlaue, tolerance vlaue of relative variation rate of the observed pseudo log-loglikelihood value, defualt as '1e-5'.}
  \item{verbose}{an optional logical value, whether output the information of the ICM-EM algorithm.}
  \item{seed}{an optional integer, the random seed in fitting PRECAST model.}
}
\details{
  Nothing 
}
\value{
  Return a \link{list} including all paramters' setting.
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
  model_set()
}