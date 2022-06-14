\name{AddAdjList}
\alias{AddAdjList}
\title{Add adjacency matrix list for a PRECASTObj object}
\description{
  Add adjacency matrix list for a PRECASTObj object to prepare for PRECAST model fitting.
}
\usage{
  AddAdjList(PRECASTObj, type="fixed_distance", platform="Visium", ...)
}
\arguments{
  \item{PRECASTObj}{a PRECASTObj object created by \link{CreatePRECASTObject}.}
  \item{type}{an optional string, specify which type of neighbors' definition. Here we provide two definition: one is "fixed_distance", the other is "fixed_number".}
  \item{platform}{a string, specify the platform of the provided data, default as "Visium". There are many platforms to be supported, including ("Visuim", "ST", "SeqFISH", 'merFISH', 'slide-seqv2', 'seqscope', "HDST"), which means there are spatial coordinates information in the metadata of PRECASTObj. The platform helps to calculate the adjacency matrix by defining the neighborhoods when type="fixed_distance" is chosen.}
  \item{...}{other arguments to be passed to \link{getAdj}, \link{getAdj_auto}  and \link{getAdj_fixedNumber} funciton.}
}
\details{
  When the type = "fixed_distance", then the spots within the Euclidean distance cutoffs from one spot are regarded as the neighbors of this spot. When the type = "fixed_number",  the K-nearest spots are regarded as the neighbors of each spot. 
}
\value{
  Return a revised PRECASTObj object by adding the adjacency matrix list.
}

\author{
  Wei Liu
}
\note{
  nothing
}
\seealso{
  \link{AddParSetting}.
}
