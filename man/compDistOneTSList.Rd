% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DistanceComputations.R
\name{compDistOneTSList}
\alias{compDistOneTSList}
\title{Compression-/Complexity-based Dissimilarity for One List}
\usage{
compDistOneTSList(tsList, symbolLimits)
}
\arguments{
\item{tsList}{A list of numeric vectors/matrixes (uni- or multi-variate
time series).}

\item{symbolLimits}{Interval boundaries which will be used to convert the
time series to a SAX representation. Should be a monotonically increasing
vector starting with -Inf and ending with +Inf.}
}
\value{
The (symmetric) dissimilarity matrix with each entry being a numeric
from the range [0,1].
}
\description{
Version of \code{\link{compDistTSList}} which computes the compression-/
complexity-based dissimilarity from each time series of a list to each time
series of the same list.
}
