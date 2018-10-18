#' Symbolic Aggregate Aproximation Limits (Standard)
#'
#' Calculates the interval boundaries for a SAX representation based on the normal
#' distribution. Same as \code{TSclust::SAX.breakpoints.table()} (not exported
#' function). This is the approach of the original authors.
#'
#' @section References:
#'
#' Lin, J., Keogh, E., Lonardi, S. & Chiu, B. (2003). A symbolic representation
#' of time series, with implications for streaming algorithms. In \emph{Proceedings
#' of the 8th acm sigmod workshop on research issues in data mining and knowledge
#' discovery} (pp. 2-11). ACM.
#'
#' @importFrom stats qnorm
#' @param alphabetSize Number of symbols (intervals).
#' @return The interval boundaries as monotonically increasing vector from negative
#' infinity to positive infinity.
#' @family SAX functions
#' @export
SAXLimitsOriginal <- function(alphabetSize = 8) {
  return(qnorm((0:alphabetSize)/alphabetSize))
}

#' Symbolic Aggregate Aproximation Limits (Flexible)
#'
#' Calculates the interval boundaries for a SAX representation based on the real
#' distribution of values in a time series list. Uses evenly-spaced quantiles on
#' the real data instead of the normal distribution proposed by the original authors.
#'
#' @section References:
#'
#' Lin, J., Keogh, E., Lonardi, S. & Chiu, B. (2003). A symbolic representation
#' of time series, with implications for streaming algorithms. In \emph{Proceedings
#' of the 8th acm sigmod workshop on research issues in data mining and knowledge
#' discovery} (pp. 2-11). ACM.
#'
#' @importFrom stats qnorm
#' @importFrom stats quantile
#' @param tsList A list of numeric vectors/matrixes (uni- or multi-variate time
#' series).
#' @param alphabetSize Number of symbols (intervals).
#' @return The interval boundaries as monotonically increasing vector from negative
#' infinity to positive infinity.
#' @family SAX functions
#' @export
SAXLimitsDataAdaptive <- function(tsList, alphabetSize = 8) {
  dataVector <- unlist(tsList)
  result <- quantile(x = dataVector, probs = (0:alphabetSize)/alphabetSize, names = FALSE)
  result[1] <- -Inf
  result[alphabetSize + 1] <- Inf
  return(result)
}

