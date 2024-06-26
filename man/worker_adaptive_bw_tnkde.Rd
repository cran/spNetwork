% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adaptive_simultaneous_tnkde_bw.R
\name{worker_adaptive_bw_tnkde}
\alias{worker_adaptive_bw_tnkde}
\title{Worker function for adaptive bandwidth for TNDE}
\usage{
worker_adaptive_bw_tnkde(
  lines,
  quad_events,
  events_loc,
  events,
  w,
  kernel_name,
  bw_net,
  bw_time,
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
\item{lines}{A feature collection of linestrings representing the underlying
network}

\item{quad_events}{a feature collection of points indicating for which events
the densities must be calculated}

\item{events_loc}{A feature collection of points representing the location of
the events}

\item{events}{A feature collection of points representing the events.
Multiple events can share the same location. They are linked by the goid
column}

\item{w}{A numeric vector with the weight of the events}

\item{kernel_name}{The name of the kernel to use (string)}

\item{bw_net}{The fixed kernel bandwidth for the network dimension. Can also be a vector
if several bandwidth must be used.}

\item{bw_time}{The fixed kernel bandwidth for the time dimension. Can also be a vector
if several bandwidth must be used.}

\item{method}{The type of NKDE to use (string)}

\item{div}{The divisor to use for the kernel. Must be "n" (the number of
events within the radius around each sampling point), "bw" (the bandwidth)
"none" (the simple sum).}

\item{digits}{The number of digits to retain from the spatial coordinates. It
ensures that topology is good when building the network. Default is 3. Too
high a precision (high number of digits) might break some connections}

\item{tol}{A float indicating the minimum distance between the events and the
lines' extremities when adding the point to the network. When points are
closer, they are added at the extremity of the lines.}

\item{sparse}{A Boolean indicating if sparse or regular matrices should be
used by the Rcpp functions. These matrices are used to store edge indices
between two nodes in a graph. Regular matrices are faster, but require more
memory, in particular with multiprocessing. Sparse matrices are slower (a
bit), but require much less memory.}

\item{max_depth}{An integer, the maximum depth to reach for continuous and
discontinuous NKDE}

\item{verbose}{A Boolean, indicating if the function should print messages
about the process.}
}
\value{
A vector with the local bandwidths or an array if bw_net and bw_time are vectors
}
\description{
The worker function to calculate Adaptive bandwidths according
to Abramson’s smoothing regimen for TNKDE with a space-time interaction (INTERNAL).
}
\examples{
#This is an internal function, no example provided
}
