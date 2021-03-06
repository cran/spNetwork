% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graph_functions.R
\name{graph_checking}
\alias{graph_checking}
\title{Topological error}
\usage{
graph_checking(lines, digits, tol)
}
\arguments{
\item{lines}{A SpatialLinesDataFrame representing the network}

\item{digits}{An integer indicating the number of digits to retain for
coordinates}

\item{tol}{A float indicating the tolerance distance to identify a dangle
node}
}
\value{
A list with two elements. The first is a SpatialPointsDataFrame
  indicating for each node of the network to which component it belongs. The
  second is a SpatialPointsDataFrame with the dangle nodes of the network.
}
\description{
A utility function to find topological errors in a network.
}
\examples{
\dontrun{
networkgpkg <- system.file("extdata", "networks.gpkg", package = "spNetwork", mustWork = TRUE)
mtl_network <- rgdal::readOGR(networkgpkg,layer="mtl_network", verbose=FALSE)
topo_errors <- graph_checking(mtl_network, 2, 2)
}
}
