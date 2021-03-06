% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{tempCor_fast}
\alias{tempCor_fast}
\title{Temporal Correlation}
\usage{
tempCor_fast(x, y)
}
\arguments{
\item{x}{1st numeric vector/time series.}

\item{y}{2nd numeric vector/time series.}
}
\value{
The temporal correlation as double from the range [-1,1].
}
\description{
Calculates the temporal correlation (correlation between the difference vectors)
as defined by Chouakria and Nagabhushan (2007). We additionally set the
correlation between two constant series to 1 and between a constant and a
non-constant one to 0 to guarantee that the result is always a proper (aka
not NaN, not infinity) number.
}
\section{References}{


Chouakria, A. D. & Nagabhushan, P. N. (2007). Adaptive dissimilarity index
for measuring time series proximity. \emph{Advances in Data Analysis and
Classification, 1}(1), 5-21.
}

