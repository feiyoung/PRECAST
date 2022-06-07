\name{AddAdjList}
\alias{AddAdjList}
\title{Add adjacency matrix list for a PRECASTObj object}
\description{
  Add adjacency matrix list for a PRECASTObj object to prepare for PRECAST model fitting.
}
\usage{
  AddAdjList(PRECASTObj, type="fixed_distance", platform="Visium", ...)
}
\arguments{
  \item{PRECASTObj}{a PRECASTObj object created by \link{CreatePRECASTObject}.}
  \item{type}{}
  \item{platform}{a string, specify the platform of the provided data, default as "Visium". There are many platforms to be supported, including ("Visuim", "ST", "SeqFISH", 'merFISH', 'slide-seqv2', 'seqscope', "HDST"), which means there are spatial coordinates information in the metadata of PRECASTObj. The platform helps to calculate the adjacency matrix by defining the neighborhoods.}
  \item{...}{other arguments to be passed to \link{getAdj}  and \link{getAdj_fixedNumber} funciton.}
}
\details{
  Nothing
}
\value{
  Return a revised PRECASTObj object by adding the adjacency matrix list.
}

\author{
  Wei Liu
}
\note{
  nothing
}
\seealso{
  \link{AddParSetting}.
}
\examples{

}