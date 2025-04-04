% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{kfunc_cpp2}
\alias{kfunc_cpp2}
\title{c++ k function 2}
\usage{
kfunc_cpp2(dist_mat, start, end, step, Lt, n, wc, wr, cross = FALSE)
}
\arguments{
\item{dist_mat}{A square matrix with the distances between points}

\item{start}{A float, the start value for evaluating the k-function}

\item{end}{A float, the last value for evaluating the k-function}

\item{step}{A float, the jump between two evaluations of the k-function}

\item{Lt}{The total length of the network}

\item{n}{The number of points}

\item{wc}{The weight of the points represented by the columns (destinations)}

\item{wr}{The weight of the points represented by the rows (origins)}

\item{cross}{A boolean indicating if we are calculating a cross k function or not (default is FALSE)}
}
\value{
A numeric vector with the values of the k function evaluated at the required distances
}
\description{
c++ k function (INTERNAL)
}
\keyword{internal}
