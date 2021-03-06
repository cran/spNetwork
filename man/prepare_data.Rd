% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nkde_execution_functions.R
\name{prepare_data}
\alias{prepare_data}
\title{Prior data preparation}
\usage{
prepare_data(samples, lines, events, w, digits, tol, agg)
}
\arguments{
\item{samples}{A spatialPointsDataFrame of the samples points}

\item{lines}{A SpatialLinesDataFrame representing the network}

\item{events}{A spatialPointsDataFrame of the events points}

\item{w}{A numeric vector representing the weight of the events}

\item{digits}{The number of digits to keep}

\item{tol}{A float indicating the spatial tolerance when snapping events on
lines}

\item{agg}{A double indicating if the points must be aggregated within a distance.
if NULL, then the points are aggregated by rounding the coordinates.}
}
\value{
the data prepared for the rest of the operations
}
\description{
A simple function to prepare data before the NKDE calculation.
}
\examples{
#This is an internal function, no example provided
}
\keyword{internal}
