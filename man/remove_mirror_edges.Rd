% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geometrical_functions_sf.R
\name{remove_mirror_edges}
\alias{remove_mirror_edges}
\title{Remove mirror edges}
\usage{
remove_mirror_edges(lines, keep_shortest = TRUE, digits = 3, verbose = TRUE)
}
\arguments{
\item{lines}{A feature collection of linestrings}

\item{keep_shortest}{A boolean, if TRUE, then the shortest line is keeped if
several lines have the same starting point and ending point. if FALSE, then the
longest line is keeped.}

\item{digits}{An integer indicating the number of digits to keep in coordinates}
}
\value{
A feature collection of linestrings with the mirror edges removed
}
\description{
Keep unique edges based on start and end point
}
\examples{
#This is an internal function, no example provided
}
\keyword{internal}
