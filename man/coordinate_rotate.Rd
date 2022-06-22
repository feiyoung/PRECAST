\name{coordinate_rotate}
\alias{coordinate_rotate}
\title{Coordinates rotation for visualization}
\description{
  Coordinates rotation for visualization. 
}
\usage{
  coordinate_rotate(pos, theta=0)
}
\arguments{
  \item{pos}{a matrix, the n-by-d coordinates, where n is the number of coordinates, d is the dimension of coordinates.}
  \item{theta}{a real number, the angle for counter-clock-wise rotation.}
}
\details{
  Nothing
}
\value{
  Return a rotated coordinate matrix.
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
    x <- 1:100
    pos <- cbind(x, sin(pi/4*x))
    oldpar <- par(mfrow = c(1,2))
    plot(pos)
    plot(coordinate_rotate(pos, 40))
    par(oldpar)
    
}