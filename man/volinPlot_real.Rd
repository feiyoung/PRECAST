\name{volinPlot}
\alias{volinPlot}
\title{Volin/boxplot plot}
\description{
  Plot volin/boxplot.
}
\usage{
  volinPlot(mat, ylabel='ARI', cols=NULL)
}
\arguments{
   \item{mat}{a matrix with columns.}
  \item{ylabel}{an optional string, the name of ylabel.}
  \item{cols}{colors used in the plot}
}
\details{
  Nothing
}
\value{
  Return a ggplot2 object.
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
   volinPlot(mat)
}