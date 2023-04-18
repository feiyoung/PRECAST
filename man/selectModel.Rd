\name{SelectModel}
\alias{SelectModel}
\alias{SelectModel.SeqK_PRECAST_Object}
\alias{SelectModel.PRECASTObj}
\title{Select best PRECAST model from candidated models}
\description{
  Select best PRECAST model from candidated models with different number of clusters.
}
\usage{
  \method{SelectModel}{SeqK_PRECAST_Object}(obj, criteria = 'MBIC',pen_const=1, return_para_est=FALSE)
  \method{SelectModel}{PRECASTObj}(obj, criteria = 'MBIC',pen_const=1, return_para_est=FALSE)
}
\arguments{
  \item{obj}{a SeqK_PRECAST_Object or  PRECASTObj object after PRECAST model fitting.}
  \item{criteria}{a string, specify the criteria used for selecting the number of clusters, supporting "MBIC", "BIC" and "AIC".}
  \item{pen_const}{an optional positive value, the adjusted constant used in the MBIC criteria. }
  \item{return_para_est}{an optional logical value, whether return the other paramters' estimators in PRECAST.}
}
\details{
  Nothing
}
\value{
  Return a revised PRECASTObj object.
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

  data(PRECASTObj)
  PRECASTObj <- SelectModel(PRECASTObj)
  
}