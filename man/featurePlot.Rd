\name{featurePlot}
\alias{featurePlot}
\title{Spatial expression heatmap}
\description{
  Plot spatial heatmap for a feature of Seurat object with spatial transcriptomics data.
}
\usage{
  featurePlot(seu, feature=NULL, cols=NULL, pt_size=1, title_size =16, quant=0.5, 
  assay='RNA' , reduction="position")
}
\arguments{
  \item{seu}{an object named "Seurat". The object of class "Seurat" must include slot "scale.data".}
  \item{feature}{an optional string, specify the name of feature to be plotted. If it is null, the first feature will be plotted.}
  \item{cols}{colors used in the plot}
  \item{pt_size}{the size of point in the spatial heatmap plot.}
  \item{title_size}{the title size used for the plot.}
  \item{quant}{the quantile value to generate the gradient color map.}
  \item{assay}{the assay selected for plot.}
  \item{reduction}{the Reduc object for plot.}
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

  library(Seurat)
  data(PRECASTObj)
  PRECASTObj <- SelectModel(PRECASTObj)
  seuInt <- IntegrateSpaData(PRECASTObj, species='unknown')
  seuInt <- ScaleData(seuInt)
  featurePlot(seuInt, assay='PRE_CAST')
  
}