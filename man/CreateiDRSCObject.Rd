\name{CreatePRECASTObject}
\alias{CreatePRECASTObject}
\title{Create the PRECAST object with preprocessing step.}
\description{
Create the PRECAST object with preprocessing step.
}
\usage{
CreatePRECASTObject(seuList,  project = "PRECAST",  gene.number=2000, 
                    selectGenesMethod='SPARK-X',numCores_sparkx=1,  
                    customGenelist=NULL, premin.spots = 20,  
                    premin.features=20, postmin.spots=15, postmin.features=15,
                              rawData.preserve=FALSE,verbose=TRUE)
}
\arguments{
\item{seuList}{a list consisting of Seurat objects, where each object is a SRT data batch. The default assay of each Seurat object will be used for data preprocessing and followed model fitting. The specified format about seuList argument can be referred to the details and example.}
\item{project}{An optional string, name of the project, default as "PRECAST".}
\item{gene.number}{an optional integer, the number of top spatially variable genes (SVGs) or highly variable  genes (HVGs) to be chosen.}
\item{selectGenesMethod}{an optional integer, the method to select genes for each sample. It supports 'SPARK-X' and 'HVGs' to select genes now. Users can provide self-selected genes using customGenelist argument.}
\item{numCores_sparkx}{an optional integer, specify the number of CPU cores in SPARK package to use when selecting spatial genes.}
\item{customGenelist}{an optional string vector, the list of user specified genes to be used for PRECAST model fitting. If this argument is given, SVGs/HVGs will not be selected.}
\item{premin.spots}{An optional integer, the features (genes) are retained in raw data filtering step with at least premin.spots number of spots, default is 20.}
\item{premin.features}{An optional integer, the locations are retained in raw data filtering step with at least premin.features number of  nonzero-count features (genes), default is 20.}
\item{postmin.spots}{An optional integer, the features (genes) are retained in filtering step after common genes selected among all data batches with at least premin.spots number of spots, default is 15.}
\item{postmin.features}{An optional integer, the locations are retained in filtering step after common genes selected among all data batches  with at least premin.features number of  nonzero-count features (genes), default is 15.}
\item{rawData.preserve}{An optional logical value, whether preserve the raw seuList data.}
\item{verbose}{whether display the message in the creating process.}
}
\details{
  seuList is a \link{list}  with Seurat object as component, and each Seruat object includes the raw expression count matrix, spatial coordinates and meta data for each data batch, where the spatial coordinates information must be saved in the metadata of Seurat, named "row" and "col" for each data batch.
}
\value{
   Returns PRECAST object prepared for PRECAST model fitting. See \link{PRECASTObj-class} for more details.
}

\examples{
  data(PRECASTObj)
  library(Seurat)
  seuList <- PRECASTObj@seulist
  ## Check the input of seuList for create PRECAST object.
  ## Check the default assay for each data batch
  lapply(seuList, DefaultAssay)
  ## Check the spatial coordinates in the meta data named "row" and "col".
  head(seuList[[1]]@meta.data)
  ## Then create PRECAST object using this seuList.
  ## For convenience, we show the  user-specified genes' list for creating PRECAST object.
  ## Users can use SVGs from SPARK-X or HVGs.
  PRECASTObj2 <- CreatePRECASTObject(seuList, 
   customGenelist= row.names(seuList[[1]]), verbose=FALSE)


}