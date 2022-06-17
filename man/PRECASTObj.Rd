\docType{class}
\name{PRECASTObj-class}
\alias{PRECASTObj-class}
\title{Each PRECASTObj object has a number of slots which store information. }
\description{
Each PRECASTObj object has a number of slots which store information. Key slots to access
are listed below.
}
\section{Slots}{

\describe{
\item{\code{seuList}}{A list with Seurat object as component, representing the raw expression count matrix, spatial coordinates and meta data for eahc data batch, where the spatial coordinates information is saved in the metadata of Seurat, named "row" and "col" for eahc data batch.}

\item{\code{seulist}}{A Seurat list after the preprocessing step in preparation for PRECAST model.}

\item{\code{AdjList}}{The adjacency matrix list for a PRECASTObj object.}

\item{\code{parameterList}}{The model parameter settings  for a PRECASTObj object}

\item{\code{resList}}{The results after fitting PRECAST models.}

\item{\code{project}}{Name of the project.}

}}