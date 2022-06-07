\name{getAdj_reg}
\alias{getAdj_reg}
\title{Calculate adjacency matrix for regular spatial coordinates.}
\description{
 Calculate adjacency matrix for regular spatial coordinates from ST or Visium platform.
}
\usage{
getAdj_reg(pos, platform= "Visium")
}
\arguments{
\item{pos}{is a n-by-d matrix of position, where n is the number of spots, and d is the dimension of coordinates.}
 \item{platform}{a string, specify the platform of the provided data, default as "Visium", and only support "ST" and "Visium" platform. }
}
\value{
 A sparse matrix containing the neighbourhood.
}
\seealso{
  \link{getAdj_auto}, \link{getAdj}, \link{getAdj_fixedNumber}.
}
