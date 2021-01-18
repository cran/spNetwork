% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geometrical_functions.R
\name{remove_loop_lines}
\alias{remove_loop_lines}
\title{Remove loops}
\usage{
remove_loop_lines(lines, digits)
}
\arguments{
\item{lines}{A SpatialLinesDataFrame}

\item{digits}{An integer indicating the number of digits to keep for the
spatial coordinates}
}
\value{
A SpatialLinesDataFrame
}
\description{
Remove from a SpatialLinesDataFrame the lines that have the
same starting and ending point.
}
\examples{
#This is an internal function, no example provided
}
\keyword{internal}