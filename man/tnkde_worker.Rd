% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/temporal_nkde_sf.R
\name{tnkde_worker}
\alias{tnkde_worker}
\title{TNKDE worker}
\usage{
tnkde_worker(
  lines,
  events_loc,
  events,
  samples_loc,
  samples_time,
  kernel_name,
  bw_net,
  bw_time,
  bws_net,
  bws_time,
  method,
  div,
  digits,
  tol,
  sparse,
  max_depth,
  verbose = FALSE
)
}
\arguments{
\item{lines}{A feature collection of linestrings with the sampling points. The
geometries must be simple Linestrings (may crash if some geometries
are invalid)}

\item{events_loc}{A feature collection of points representing the aggergated events on the
network. The points will be snapped on the network.}

\item{events}{A feature collection of points representing the base events on the
network}

\item{samples_loc}{A feature collection of points representing the locations for
which the densities will be estimated.}

\item{samples_time}{A numeric vector representing when each density will be
estimated}

\item{kernel_name}{The name of the kernel to use}

\item{bw_net}{The global network kernel bandwidth}

\item{bw_time}{The global time kernel bandwidth}

\item{bws_net}{The network kernel bandwidth (in meters) for each event}

\item{bws_time}{The time bandwidth for each event}

\item{method}{The method to use when calculating the NKDE, must be one of
simple / discontinuous / continuous (see details for more information)}

\item{div}{The divisor to use for the kernel. Must be "n" (the number of
events within the radius around each sampling point), "bw" (the bandwidth)
"none" (the simple sum).}

\item{digits}{The number of digits to keep in the spatial coordinates. It
ensures that topology is good when building the network. Default is 3}

\item{tol}{When adding the events and the sampling points to the network,
the minimum distance between these points and the lines extremities. When
points are closer, they are added at the extremity of the lines.}

\item{sparse}{A Boolean indicating if sparse or regular matrices should be
used by the Rcpp functions. Regular matrices are faster, but require more
memory and could lead to error, in particular with multiprocessing. Sparse
matrices are slower, but require much less memory.}

\item{max_depth}{When using the continuous and discontinuous methods, the
calculation time and memory use can go wild  if the network has a lot of
small edges (area with a lot of intersections and a lot of events). To
avoid it, it is possible to set here a maximum depth. Considering that the
kernel is divided at intersections, a value of 8 should yield good
estimates. A larger value can be used without problem for the discontinuous
method. For the continuous method, a larger value will strongly impact
calculation speed.}

\item{verbose}{A Boolean, indicating if the function should print messages
about the process.}
}
\value{
A numeric matrix with the nkde values
}
\description{
The worker function for tnkde and tnkde.mc
}
\examples{
#This is an internal function, no example provided
}
\keyword{internal}
