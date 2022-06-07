\name{boxPlot}
\alias{boxPlot}
\title{Boxplot for a matrix}
\description{
  Boxplot for a matrix.
}
\usage{
  boxPlot(mat, ylabel='ARI', cols=NULL, ...)
}
\arguments{
  \item{mat}{a matrix with columns.}
  \item{ylabel}{an optional string, the name of ylabel.}
  \item{cols}{colors used in the plot}
  \item{...}{Other parameters passed to geom_boxplot.}
}
\details{
  Nothing
}
\value{
  Return a ggplot2 object.
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
   mat <- matrix(runif(100*3, 0.6, 1), 100, 3)
   colnames(mat) <- paste0("Method", 1:3)
   boxPlot(mat)
}