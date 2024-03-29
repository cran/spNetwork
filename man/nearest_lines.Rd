% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geometrical_functions_sf.R
\name{nearest_lines}
\alias{nearest_lines}
\title{Nearest line for points}
\usage{
nearest_lines(points, lines, snap_dist = 300, max_iter = 10)
}
\arguments{
\item{points}{A feature collection of points}

\item{lines}{A feature collection of linestrings}

\item{snap_dist}{A distance (float) given to find for each point its
nearest line in a spatial index. A too big value will produce
unnecessary distance calculations and a too short value will lead to
more iterations to find neighbours. In extrem cases, a too short value
could lead to points not associated with lines (index = -1).}

\item{max_iter}{An integer indicating how many iteration the search
algorithm must perform in the spatial index to find lines close to a
point. At each iteration, the snap_dist is doubled to find candidates.}
}
\description{
Find for each point its nearest LineString
}
\examples{
# this is an internal function, no example provided
}
\keyword{internal}
