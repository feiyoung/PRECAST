\name{CreatePRECASTObject}
\alias{CreatePRECASTObject}
\title{Create the PRECAST object with preprocessing step.}
\description{
Create the PRECAST object with preprocessing step.
}
\usage{
CreatePRECASTObject(seuList,  project = "PRECAST", numCores_sparkx=1, 
                gene.number=2000,customGenelist=NULL, premin.spots = 20,  
                premin.features=20, postmin.spots=15, postmin.features=15,
                rawData.preserve=FALSE,verbose=TRUE)
}
\arguments{
\item{seuList}{a list consisting of Seurat objects, where each object is a SRT data batch.}
\item{project}{An optional string, name of the project, default as "PRECAST".}
\item{numCores_sparkx}{an optional integer, specify the number of CPU cores in SPARK package to use when selecting spatial genes.}
\item{gene.number}{an optional integer, the number of top spatially variable genes (SVGs) to be chosen.}
\item{customGenelist}{an optional string vector, the list of user specified genes to be used for PRECAST model fitting. If this argument is given, SVGs will not be selected.}
\item{premin.spots}{An optional integer, the features (genes) are retained in raw data filtering step with at least premin.spots number of spots, default is 20.}
\item{premin.features}{An optional integer, the locations are retained in raw data filtering step with at least premin.features number of  nonzero-count features (genes), default is 20.}

\item{postmin.spots}{An optional integer, the features (genes) are retained in filtering step after common genes selected among all data batches with at least premin.spots number of spots, default is 15.}
\item{postmin.features}{An optional integer, the locations are retained in filtering step after common genes selected among all data batches  with at least premin.features number of  nonzero-count features (genes), default is 15.}
\item{rawData.preserve}{An optional logical value, whether preserve the raw seuList data.}
\item{verbose}{whether display the message in the creating process.}
}
\value{
   Returns PRECAST object prepared for PRECAST model fitting.
}

\examples{
  data(PRECASTObj)
  seuList <- PRECASTObj@seulist
  PRECASTObj2 <- CreatePRECASTObject(seuList, 
   customGenelist= row.names(seuList[[1]]), verbose=FALSE)


}