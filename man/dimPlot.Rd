\name{dimPlot}
\alias{dimPlot}
\title{Low-dimensional embeddings' plot}
\description{
  Low-dimensional embeddings' plot colored by a specified meta data in the Seurat object.
}
\usage{
  dimPlot(seuInt, item=NULL, reduction=NULL, point_size=1,text_size=16, 
                    cols=NULL,font_family='', border_col="gray10",
                    fill_col="white", ...)
}
\arguments{
  \item{seuInt}{an object named "Seurat". }
  \item{item}{the item used for coloring the plot in the meta data of seuInt object.}
  \item{reduction}{the reduction used for plot in the seuInt object. If reduction is null, the last added one is used for plotting.}
  \item{point_size}{the size of point in the scatter plot.}
  \item{text_size}{the text size in the plot.}
  \item{cols}{colors used in the plot}
  \item{font_family}{the font family used for the plot.}
  \item{border_col}{the border color in the plot.}
  \item{fill_col}{the  color used in backgroup.}
  \item{...}{other arguments passed to \code{\link{plot_scatter}}}.
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
  data(PRECASTObj)
  PRECASTObj <- selectModel.PRECASTObj(PRECASTObj)
  seuInt <- IntegrateSpaData(PRECASTObj, species='unknown')
  dimPlot(seuInt, reduction = 'PRECAST')
  ## or use the Seurat::DimPlot(seuInt, reduction = 'PRECAST')
  
}