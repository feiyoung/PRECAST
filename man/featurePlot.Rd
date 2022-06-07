\name{featurePlot}
\alias{featurePlot}
\title{Spatial heatmap}
\description{
  Plot spatial heatmap for a Seurat object with spatial transcriptomics data.
}
\usage{
  featurePlot(seu, feature=NULL, cols=NULL, pt_size=1, title_size =16, quant=0.5, 
  assay='RNA' , reduction="position")
}
\arguments{
  \item{seu}{an object named "Seurat". The object of class "Seurat" must include slot "scale.data".}
  \item{feature}{an optional positive integer, specify the number of features to be extracted.}
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
  For Seurat object, return a Seurat object. For objcet "matrix" and "dgCMatrix", return a object "matrix" with rownames same as the colnames of \code{X}, and colnames "WPCA1" to "WPCAq".
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
\dontrun{
  data(PRECASTObj)
  PRECASTObj <- selectModel(PRECASTObj)
  seuInt <- IntegrateSpaData(PRECASTObj, species='unknown')
  seuInt <- ScaleData(seuInt)
  featurePlot(seuInt, assay='PRE_CAST')
  }
}