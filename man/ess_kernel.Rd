% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernel_functions_sf.R
\name{ess_kernel}
\alias{ess_kernel}
\title{Worker for simple NKDE algorithm}
\usage{
ess_kernel(graph, y, bw, kernel_func, ok_samples, nodes, ok_edges, N)
}
\arguments{
\item{graph}{a graph object from igraph representing the network}

\item{y}{the index of the actual event}

\item{bw}{a float indicating the kernel bandwidth (in meters)}

\item{kernel_func}{a function obtained with the function select_kernel}

\item{ok_samples}{a a feature collection of points representing the sampling points. The
samples must be snapped on the network. A column edge_id must indicate for
each sample on which edge it is snapped.}

\item{nodes}{a a feature collection of points representing the nodes of the network}

\item{ok_edges}{a a feature collection of linestrings representing the edges of the network}
}
\description{
The worker function to perform the simple nkde.
}
\examples{
#This is an internal function, no example provided
}
\keyword{internal}
