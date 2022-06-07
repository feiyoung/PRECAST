\name{ridgePlot}
\alias{ridgePlot}
\title{ridgePlot for features in Seurat}
\description{
  ridgePlot for features in Seurat object.
}
\usage{
  ridgePlot(object, features, group.names='Domain' , 
     cols = NULL, idents = NULL, sort = FALSE, 
     assay = NULL, group.by = NULL, y.max = NULL, same.y.lims = FALSE, 
     log = FALSE, ncol = NULL, slot = "data", stack = FALSE, 
     combine = TRUE, fill.by = "feature")
}
\arguments{
  \item{object}{an object named "Seurat".}
  \item{features}{an optional positive integer, specify the number of features to be extracted.}
  \item{group.names}{the size of point in the scatter plot.}
  \item{cols}{colors to use for plotting.}
  \item{idents}{which classes to include in the plot (default is all)}
  \item{sort}{sort identity classes (on the x-axis) by the average expression of the attribute being potted, can also pass 'increasing' or 'decreasing' to change sort direction.}
  \item{assay}{name of assay to use, defaults to the active assay.}
  \item{group.by }{group (color) cells in different ways.}
  \item{ y.max }{maximum y axis value}
  \item{same.y.lims}{set all the y-axis limits to the same values}
  \item{log}{plot the feature axis on log scale}
  \item{ncol}{number of columns if multiple plots are displayed}
  \item{slot}{use non-normalized counts data for plotting}
  \item{stack}{horizontally stack plots for each feature}
  \item{combine}{combine plots into a single patchworked ggplot object. If FALSE, return a list of ggplot}
  \item{fill.by}{color violins/ridges based on either 'feature' or 'ident'}
}
\details{
  Nothing
}
\value{
  Return a single patchworked ggplot object. If coombine = FALSE, return a list of ggplot.
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
  ridgePlot(seuInt, features=row.names(seuInt)[1:4], ncol=2)
  }
}