\name{doHeatmap}
\alias{doHeatmap}
\title{Heatmap for spots-by-feature matrix}
\description{
  Plot  heatmap for a Seurat object with expressioin data.
}
\usage{
  doHeatmap(seu, features=NULL, cell_label='Cell type', grp_label = FALSE,
                      pt_size=4, grp_color=NULL, ...)
}
\arguments{
  \item{seu}{an object named "Seurat". The object of class "Seurat" must include slot "scale.data".}
  \item{features}{an optional string vector, the features to be plotted.}
  \item{cell_label}{an optional string, the name of legend.}
  \item{grp_label}{an optional logical value, whether display the group names.}
  \item{pt_size}{the point size used in the plot}
  \item{grp_color}{the colors to use for the group color bar.}
  \item{...}{Other paramters passed to \link{DoHeatmap}.}
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
  \link{featurePlot}
}
\examples{

  library(Seurat)
  data(PRECASTObj)
  PRECASTObj <- selectModel(PRECASTObj)
  seuInt <- IntegrateSpaData(PRECASTObj, species='unknown')
  seuInt <- ScaleData(seuInt)
  doHeatmap(seuInt, features=row.names(seuInt)[1:5])
  
}