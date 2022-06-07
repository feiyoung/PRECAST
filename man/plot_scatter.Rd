\name{plot_scatter}
\alias{plot_scatter}
\title{Spatial heatmap}
\description{
  Plot spatial heatmap for a Seurat object with spatial transcriptomics data.
}
\usage{
  plot_scatter(embed_use, meta_data, label_name, xy_names=c('tSNE1', 'tSNE2'), no_guides = FALSE, 
    palette_use = tableau_color_pal()(10), 
    point_size = 0.5, point_alpha=1, 
    base_size = 12, do_points = TRUE, do_density = FALSE, border_col='gray',
    legend_pos='right', legend_dir='vertical')
}
\arguments{
  \item{embed_use}{an object named "Seurat", "maxtrix" or "dgCMatrix". The object of class "Seurat" must include slot "scale.data".}
  \item{meta_data}{an optional positive integer, specify the number of features to be extracted.}
  \item{label_name}{the size of point in the scatter plot.}
  \item{xy_names}{the text size in the plot.}
  \item{no_guides}{colors used in the plot}
  \item{palette_use}{the font family used for the plot.}
  \item{point_size}{the point size of scatter plot.}
  \item{point_alpha}{the transparency of the plot.}
  \item{base_size}{the base text size.}
  \item{do_points}{Plot point.}
  \item{do_density}{Plot density plot}
  \item{border_col}{the border color in the plot.}
  \item{legend_pos}{the position of legend.}
  \item{legend_dir}{the direction of legend.}
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
  embed_use <- cbind(1:100, sin((1:100)*pi/2))
  meta_data <- data.frame(cluster=factor(rep(1:2, each=50)))
  plot_scatter(embed_use, meta_data, label_name='cluster')
}