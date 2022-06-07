\name{AddParSetting}
\alias{AddParSetting}
\title{Add model settings  for a PRECASTObj object}
\description{
  The main interface function provides serveral PRECAST submodels, so a model setting is required to specified in advance  for a PRECASTObj object.
}
\usage{
  AddParSetting(PRECASTObj, ...)
}
\arguments{
  \item{PRECASTObj}{a PRECASTObj object created by \link{CreatePRECASTObject}.}
  \item{...}{other arguments to be passed to \link{model_set} funciton.}
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
  PRECASTObj <-AddParSetting(PRECASTObj)
  PRECASTObj@parameterList
}