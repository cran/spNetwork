% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nkde_execution_functions.R
\name{aggregate_points}
\alias{aggregate_points}
\title{Events aggregation}
\usage{
aggregate_points(points, maxdist)
}
\arguments{
\item{points}{The SpatialPointsDataFrame to contract (must have a weight column)}

\item{maxdist}{The distance to use}
}
\value{
A new SpatialPointsDataFrame
}
\description{
Function to aggregate points within a radius.
}
\examples{
#This is an internal function, no example provided
}
\keyword{internal}
