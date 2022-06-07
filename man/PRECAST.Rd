\name{PRECAST}
\alias{PRECAST}
\title{Fit a PRECAST model}
\description{
  Fit a PRECAST model.
}
\usage{
  PRECAST(PRECASTObj, K=NULL, q= 15)
}
\arguments{
  \item{PRECASTObj}{an object named "PRECASTObj". The object PRECASTObj is created by \link{CreatePRECASTObject}.}
  \item{K}{An optional integer or integer vector, specify the candidates of number of clusters. if \code{K=NULL}, it will be set to 4~12.}
  \item{q}{An optional integer, specify the number of low-dimensional embeddings to extract in  PRECAST.}
}
\details{
  The model fitting results are saved in the slot of resList.
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

}