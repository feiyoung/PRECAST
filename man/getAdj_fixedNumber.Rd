\name{getAdj_fixedNumber}
\alias{getAdj_fixedNumber}
\title{Calculate adjacency matrix by user-specified number of neighbors}
\usage{
getAdj_fixedNumber(pos, number=6)
}
\arguments{
\item{pos}{is a n-by-d matrix of position, where n is the number of spots, and d is the dimension of coordinates.}

\item{number}{is the number of neighbors of each spot. Euclidean distance to decide whether a spot is an neighborhood of another spot.}
}
\value{
 A sparse matrix containing the neighbourhood.
}
\seealso{
  \link[DR.SC]{getAdj}, \link[DR.SC]{getAdj_auto}.
}
\description{
 an efficient function to find the neighborhood based on the matrix of position and a  user-specified number of neighbors of each spot.
}