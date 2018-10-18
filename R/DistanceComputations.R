#' Average Time Series
#'
#' Averages a list of equal-length uni- or multivariate time series based on
#' dissimilarities (which are normalized, transformed to similarities and used
#' as weights), puting more weight on series which have a higher similarity
#' (lower dissimilarity).
#'
#' @param tsList A list of numeric vectors/matrices (uni-/multivariate time series).
#' @param dissimilarities A numeric vector of dissimilarity values (same length
#' as the list of time series).
#' @param dissWeightConversion Strategy for converting similarities to
#' dissimilarities. Might be "minmax" (linear), "gaussian" or "laplacian".
#' @export
averageTimeSeriesList <- function(tsList, dissimilarities, dissWeightConversion = "minmax") {
  maxDiss <- max(dissimilarities)
  minDiss <- min(dissimilarities)
  if (maxDiss == minDiss) {
    tsWeights <- rep(1/length(tsList), times = length(tsList))
  } else {
    if (dissWeightConversion == "minmax") {
      similarities <- 1 - dissimilarities / maxDiss
    } else if (dissWeightConversion == "laplacian") {
      # normalization necessary or all similarities become really low or even 0
      # because of exponential function (if dissimilarities high)
      similarities <- exp(-(dissimilarities - minDiss))
    }  else if (dissWeightConversion == "gaussian") {
      similarities <- exp(-(dissimilarities - minDiss)^2)
    } else {
      stop("Unknown dissWeightConversion in averageTimeSeriesList()")
    }
    tsWeights <- similarities / sum(similarities) # relative
  }
  if (is.matrix(tsList[[1]])) {
    return(averageTimeSeriesMult_fast(tsList, tsWeights))
  } else {
    return(averageTimeSeries_fast(tsList, tsWeights))
  }
}

#' Normalized Cross-Correlation Function
#'
#' Computes the cross-correlation function between two vectors and normalizes
#' it with the l2 norms of the vectors (so values are between -1 and 1). Internally,
#' convolution and FFT are used as described by Paparrizos and Gravano (2015).
#'
#' @section References:
#' Paparrizos, J. & Gravano, L. (2015). K-shape: Efficient and accurate clustering
#' of time series. In \emph{Proceedings of the 2015 acm sigmod international
#' conference on management of data} (pp. 1855–1870). ACM.
#'
#' @importFrom stats convolve
#' @param x 1st numeric vector/time series.
#' @param y 2nd numeric vector/time series.
#' @return The normalized cross-correlation as vector (length equals
#' \emph{|x| + |y| - 1}.
#' @family cross-correlation functions
#' @export
crossCorNormalized <- function(x,y) {
  nx <- length(x)
  ny <- length(y)
  # Constant-0 series would result in NaNs during normalization
  if (all(x == 0)) {
    if (all(y == 0)) {
      # produce same result as two series of another constant
      nShort <- min(nx, ny)
      conv <- rep(nShort, nx + ny - 1)
      conv[1:(nShort - 1)] <- 1:(nShort - 1)
      conv[(length(conv) - nShort + 2):length(conv)] <- (nShort - 1):1
      return(conv / sqrt(nx * ny))
    } else {
      return(rep(0, times = nx + ny - 1))
    }
  } else if (all(y == 0)) {
    return(rep(0, times = nx + ny - 1))
  }
  # Pad time series to remove circular influence; convolve(type = "open") applied
  # to unpadded series would yield same result, but is a lot slower for some reason
  # (probably because no "clever" padding to a power of 2 in R implementation)
  nConv <- 2^ceiling(log2(2*max(nx,ny)))
  convolution <- convolve(x = c(x, rep(0, times = nConv - nx)),
                          y = c(y, rep(0, times = nConv - ny)))
  if (ny == 1) {
    relevantConvValues <- convolution[1:nx]
  } else {
    relevantConvValues <- c(convolution[(nConv - ny + 2):nConv], convolution[1:nx])
  }
  return(relevantConvValues / (l2Norm_fast(x)*l2Norm_fast(y)))
}

#' Shape-Based Distance
#'
#' Computes the dissimilarity as \emph{1 - maximum normalized cross-correlation
#' coefficient} as described by Paparrizos and Gravano (2015). Multi-variate
#' time series are handled by calculating the cross-correlation between corresponding
#' attributes in \code{x} and \code{y}, averaging over attributes and then taking
#' the average cross-correlation at the lag which maximizes it.
#'
#' @section References:
#' Paparrizos, J. & Gravano, L. (2015). K-shape: Efficient and accurate clustering
#' of time series. In \emph{Proceedings of the 2015 acm sigmod international
#' conference on management of data} (pp. 1855–1870). ACM.
#'
#' @param x 1st numeric vector/time series.
#' @param y 2nd numeric vector/time series.
#' @return The dissimilarity as numeric from the range [0,2].
#' @family cross-correlation functions
#' @export
shapeBasedDistance <- function(x,y) {
  if (length(x) == length(y) && all(x == y)) {
    return(0); # important property of a metric (numerical result slightly greater 0)
  }
  # max(0, ...) because of numerical imprecision
  if (is.matrix(x) && is.matrix(y)) { # maximize (average) cross-correlation over lags
    return(max(0, 1 - max(rowMeans(sapply(1:ncol(x), function(j) {
      crossCorNormalized(x[,j], y[,j])
    })))))
  } else {
    return(max(0, 1 - max(crossCorNormalized(x,y))))
  }
}

#' L2 Distance
#'
#' Computes the standard Euclidean distance.
#'
#' @param x 1st numeric vector/time series.
#' @param y 2nd numeric vector/time series.
#' @return The L2 norm as numeric.
#' @family L_p distances
#' @export
l2Dist <- function(x,y) {
  return(sqrt(sum((x - y)^2)))
}

#' Formula for the Complexity-Based Dissimilarity
#'
#' Method which contains the formula used in \code{\link{compDist}} and
#' \code{\link{compDistTSList}}. Considering the complexity to be a kind of
#' entropy measure, our formula is similar to the normalized Variation of
#' Information described by Meila (2003) and Wu, Xiong and Chen (2009), setting
#' the joint entropy in relation to the single entropies. This results in a value
#' from the interval [0,1], compared to (0.5,1] in the formula of Keogh et al.
#' (2007). Our measure is symmetric.
#'
#' @section References:
#'
#' Keogh, E., Lonardi, S., Ratanamahatana, C. A., Wei, L., Lee, S.-H. & Handley, J.
#' (2007). Compression-based data mining of sequential data. \emph{Data Mining
#' and Knowledge Discovery, 14}(1), 99–129.
#'
#' Li, M., Badger, J. H., Chen, X., Kwong, S., Kearney, P. & Zhang, H. (2001).
#' An information-based sequence distance and its application to whole
#' mitochondrial genome phylogeny. \emph{Bioinformatics, 17}(2), 149–154.
#'
#' Meila, M. (2003). Comparing clusterings by the variation of information. In
#' B. Schölkopf & M. K. Warmuth (Eds.), \emph{Learning theory and kernel machines:
#' 16th annual conference on learning theory and 7th kernel workshop, colt/kernel
#' 2003, washington, dc, usa, august 24-27, 2003. proceedings} (pp. 173-187).
#' Springer Berlin Heidelberg.
#'
#' Wu, J., Xiong, H. & Chen, J. (2009). Adapting the right measures for k-means
#' clustering. In \emph{Proceedings of the 15th acm sigkdd international conference
#' on knowledge discovery and data mining} (pp. 877-886). ACM.
#'
#' @param xLength Length of the first string/time series after compression.
#' @param yLength Length of the second string/time series after compression.
#' @param xyLength Length of the concatenation of first and second string/time
#' series after compression.
#' @param yxLength Length of the concatenation of second and first string/time
#' series after compression.
#' @return The dissimilarity as numeric from the range [0,1].
calcCompDist <- function(xLength, yLength, xyLength, yxLength) {
  # rescale to [0,1]; max(0, ...) needed because xyLength might be smaller
  # than xLength and yLength (e.g. try with c("1", "1", "1") and c("3", "1", "1"))
  return(min(1, max(0, 2 * min(xyLength, yxLength) / (xLength + yLength) - 1)))
}

#' Compression-/Complexity-based Dissimilarity
#'
#' Dissimilarity based on the length of the compressed single as well as concatenated
#' time series as described by Li et al. (2001). Time series are represented with
#' SAX first and then zipped, both according to Keogh et al. (2007). As an
#' improvement, the dissimilarity is scaled to the interval [0,1] (before: (0.5,1])
#' and made symmetric. Multi-variate time series are handled by attribute
#' concatenation.
#'
#' @section References:
#'
#' Keogh, E., Lonardi, S., Ratanamahatana, C. A., Wei, L., Lee, S.-H. & Handley, J.
#' (2007). Compression-based data mining of sequential data. \emph{Data Mining
#' and Knowledge Discovery, 14}(1), 99–129.
#'
#' Li, M., Badger, J. H., Chen, X., Kwong, S., Kearney, P. & Zhang, H. (2001).
#' An information-based sequence distance and its application to whole
#' mitochondrial genome phylogeny. \emph{Bioinformatics, 17}(2), 149–154.
#'
#' @param x 1st numeric vector/matrix (uni- or multi-variate time series).
#' @param y 2nd numeric vector/matrix (uni- or multi-variate time series).
#' @param symbolCount Number of SAX symbols. Boundaries for the intervals will
#' be determined based on the standard normal distribution. As an alternative,
#' you can supply the boundaries directly.
#' @param symbolLimits Interval boundaries which will be used to convert the
#' time series to a SAX representation. Should be a monotonically increasing
#' vector starting with -Inf and ending with +Inf. The parameter
#' \code{symbolCount} is ignored if you supply a value here.
#' @return The dissimilarity as numeric from the range [0,1].
#' @family compression-based distances
#' @export
compDist <- function(x,y, symbolCount = 8, symbolLimits = NULL) {
  if (is.null(symbolLimits)) {
    symbolLimits <- SAXLimitsOriginal(alphabetSize = symbolCount)
  }
  xString <- as.character(SAX_fast(as.numeric(x), symbolLimits))
  yString <- as.character(SAX_fast(as.numeric(y), symbolLimits))
  if (length(xString) == length(yString) && all(xString == yString)) {
    return(0)
  }
  xyString <- c(xString, yString)
  yxString <- c(yString, xString)
  xLength <- length(memCompress(xString, type = "gzip"))
  yLength <- length(memCompress(yString, type = "gzip"))
  xyLength <- length(memCompress(xyString, type = "gzip"))
  yxLength <- length(memCompress(yxString, type = "gzip"))
  return(calcCompDist(xLength, yLength, xyLength, yxLength))
}

#' Compression-/Complexity-based Dissimilarity for Two Lists
#'
#' Version of \code{\link{compDistTSList}} which computes the compression-/
#' complexity-based dissimilarity from each time series of the first list to
#' each time series of the second list.
#'
#' @param tsList1 A list of numeric vectors/matrixes (uni- or multi-variate
#' time series).
#' @param tsList2 A list of numeric vectors/matrixes (uni- or multi-variate
#' time series). Should have the same data structure as the first parameter.
#' @param symbolLimits Interval boundaries which will be used to convert the
#' time series to a SAX representation. Should be a monotonically increasing
#' vector starting with -Inf and ending with +Inf.
#' @return The dissimilarity matrix with each entry being a numeric from the
#' range [0,1].
compDistTwoTSLists <- function(tsList1, tsList2, symbolLimits) {
  stringRepList1 <- lapply(tsList1, function(x) as.character(SAX_fast(as.numeric(x), symbolLimits)))
  stringRepList2 <- lapply(tsList2, function(x) as.character(SAX_fast(as.numeric(x), symbolLimits)))
  stringCompLength1 <- sapply(stringRepList1, function(x) length(memCompress(x, type = "gzip")))
  stringCompLength2 <- sapply(stringRepList2, function(x) length(memCompress(x, type = "gzip")))
  rowCount <- length(tsList1)
  colCount <- length(tsList2)
  distMat <- matrix(0, nrow = rowCount, ncol = colCount)
  for (i in 1:rowCount) {
    for (j in 1:colCount) {
      if (length(stringRepList1[[i]]) == length(stringRepList2[[j]]) &&
          all(stringRepList1[[i]] == stringRepList2[[j]])) {
        distance <- 0
      } else {
        xyString <- c(stringRepList1[[i]], stringRepList2[[j]])
        yxString <- c(stringRepList2[[j]], stringRepList1[[i]])
        xyLength <- length(memCompress(xyString, type = "gzip"))
        yxLength <- length(memCompress(yxString, type = "gzip"))
        distance <- calcCompDist(stringCompLength1[i], stringCompLength2[j],
                                 xyLength, yxLength)
      }
      distMat[i,j] <- distance
    }
  }
  return(distMat)
}

#' Compression-/Complexity-based Dissimilarity for One List
#'
#' Version of \code{\link{compDistTSList}} which computes the compression-/
#' complexity-based dissimilarity from each time series of a list to each time
#' series of the same list.
#'
#' @param tsList A list of numeric vectors/matrixes (uni- or multi-variate
#' time series).
#' @param symbolLimits Interval boundaries which will be used to convert the
#' time series to a SAX representation. Should be a monotonically increasing
#' vector starting with -Inf and ending with +Inf.
#' @return The (symmetric) dissimilarity matrix with each entry being a numeric
#' from the range [0,1].
compDistOneTSList <- function(tsList, symbolLimits) {
  stringRepList <- lapply(tsList, function(x) as.character(SAX_fast(as.numeric(x), symbolLimits)))
  stringCompLength <- sapply(stringRepList, function(x) length(memCompress(x, type = "gzip")))
  tsCount <- length(tsList)
  distMat <- matrix(0, nrow = tsCount, ncol = tsCount)
  for (i in 1:(tsCount - 1)) {
    for (j in (i + 1):tsCount) {
      if (length(stringRepList[[i]]) == length(stringRepList[[j]]) &&
          all(stringRepList[[i]] == stringRepList[[j]])) {
        distance <- 0
      } else {
        xyString <- c(stringRepList[[i]], stringRepList[[j]])
        yxString <- c(stringRepList[[j]], stringRepList[[i]])
        xyLength <- length(memCompress(xyString, type = "gzip"))
        yxLength <- length(memCompress(yxString, type = "gzip"))
        distance <- calcCompDist(stringCompLength[i], stringCompLength[j],
                                 xyLength, yxLength)
      }
      distMat[i,j] <- distance
      distMat[j,i] <- distance
    }
  }
  return(distMat)
}

#' Compression-/Complexity-based Dissimilarity
#'
#' Version of \code{\link{compDist}} which operates on a list of time series and
#' returns a distance matrix instead of a single distance. Saves computation time
#' compared to naive n^2 calls of the original function by computing the SAX
#' representations and single time series compression lengths only once for each
#' time series (not in every distance computation).
#'
#' @section References:
#'
#' Keogh, E., Lonardi, S., Ratanamahatana, C. A., Wei, L., Lee, S.-H. & Handley, J.
#' (2007). Compression-based data mining of sequential data. \emph{Data Mining
#' and Knowledge Discovery, 14}(1), 99–129.
#'
#' Li, M., Badger, J. H., Chen, X., Kwong, S., Kearney, P. & Zhang, H. (2001).
#' An information-based sequence distance and its application to whole
#' mitochondrial genome phylogeny. \emph{Bioinformatics, 17}(2), 149–154.
#'
#' @param tsList 1) A list of numeric vectors/matrixes (uni- or multi-variate
#' time series). The dissimilarity of the list to itself (each time series to
#' each time series) will be computed, resulting in a symmetric dissimilarity
#' matrix. 2) A list with two components, each being a list of numeric vectors/
#' matrixes (uni- or multi-variate time series). The dissimilarity of each time
#' series from the 1st component to each time series from the 2nd component will
#' be computed.
#' @param symbolCount Number of SAX symbols. Boundaries for the intervals will
#' be determined based on the standard normal distribution. As an alternative,
#' you can supply the boundaries directly.
#' @param symbolLimits Interval boundaries which will be used to convert the
#' time series to a SAX representation. Should be a monotonically increasing
#' vector starting with -Inf and ending with +Inf. The parameter
#' \code{symbolCount} is ignored if you supply a value here.
#' @return The dissimilarity matrix with each entry being a numeric from the
#' range [0,1].
#' @family compression-based distances
#' @export
compDistTSList <- function(tsList, symbolCount = 8, symbolLimits = NULL) {
  if (length(tsList) == 2 && is.list(tsList[[1]]) && is.list(tsList[[2]])) {
    if (is.null(symbolLimits)) {
      # use only second list because of compatibility in cluster assignment
      # (2nd list are the original time series)
      symbolLimits <- SAXLimitsDataAdaptive(tsList = tsList[[2]],
                                            alphabetSize = symbolCount)
    }
    return(compDistTwoTSLists(tsList1 = tsList[[1]], tsList2 = tsList[[2]],
                              symbolLimits = symbolLimits))
  } else {
    if (is.null(symbolLimits)) {
      symbolLimits <- SAXLimitsDataAdaptive(tsList = tsList,
                                            alphabetSize = symbolCount)
    }
    return(compDistOneTSList(tsList = tsList, symbolLimits = symbolLimits))
  }
}

#' Min-Max Normalization
#'
#' Normalizes a time series to [0,1]. Multi-variate time series (matrices) are
#' normalized per column.
#'
#' @param x A numeric vector/matrix (uni- or multivariate time series).
#' @return Normalized time series.
#' @export
minMaxNormalize <- function(x) {
  if (is.matrix(x)) {
    return(apply(x, 2, minMaxNormalize))
  }
  tsMin <- min(x)
  tsMax <- max(x)
  if (tsMin == tsMax) {
    return(rep(0, times = length(x)))
  } else {
    return((x - tsMin) / (tsMax - tsMin))
  }
}

#' Permutation Distribution Entropy
#'
#' More stable version of the entropy measure which is used by Brandmaier (2015)
#' to determine the optimal embedding length in permutation distribution clustering.
#' Mainly copies the (internal) function \code{pdc::codebook.entropy}, but also
#' handles the case where only one pattern occurs in the time series.
#'
#' Encodes time series subsequences of the requested length, counts the relative
#' frequency of each ordinal pattern and computes the entropy on this distribution.
#' To obtain an unbiased estimate, Brandmaier proposes to divide by the log of
#' non-zero-bin count.
#'
#' @section References:
#'
#' Brandmaier, A. M. (2011). \emph{Permutation distribution clustering and structural
#' equation model trees} (Doctoral dissertation, Universität des Saarlandes, Saarbrücken).
#'
#' Brandmaier, A. M. (2015). pdc: An R Package for Complexity-Based Clustering of
#' Time Series. \emph{Journal of Statistical Software, 67}(5), 1-23.
#'
#' @param timeSeries A numeric vector/time series (minimum length: 2 elements).
#' @param subSequenceLength Number of elements in subsequences which are used to
#' count patterns.
#' @return The normalized entropy as double.
#' @family PDC functions
pdcEntropy <- function(timeSeries, subSequenceLength) {
  tsCodebook <- codebook(timeSeries, subSequenceLength)
  tsCodebook <- tsCodebook[tsCodebook != 0]
  if (length(tsCodebook) == 1) {
    return(0)
  } else {
    return(-sum(tsCodebook * log(tsCodebook))/log(length(tsCodebook)))
  }
}

#' Permutation Distribution Entropy Heuristic
#'
#' Simplified and stabilized version of the entropy heuristic used by Brandmaier (2015)
#' (\code{pdc::entropyHeuristic}) to determine the subsequence length in permutation
#' distribution clustering. Shortens the code a bit and can handle time series
#' with only one permutation pattern.
#'
#' @section References:
#'
#' Brandmaier, A. M. (2011). \emph{Permutation distribution clustering and structural
#' equation model trees} (Doctoral dissertation, Universität des Saarlandes, Saarbrücken).
#'
#' Brandmaier, A. M. (2015). pdc: An R Package for Complexity-Based Clustering of
#' Time Series. \emph{Journal of Statistical Software, 67}(5), 1-23.
#'
#' @param tsList A list of numeric vectors/time series (which should each have a
#' minimum length of two).
#' @return The optimal subsequence/pattern length as integer from the range [2,7].
#' @family PDC functions
#' @export
pdcEntropyHeuristic <- function(tsList) {
  minTSLength <- min(sapply(tsList, length))
  if (minTSLength == 1) {
    stop("Minimum time series length of 2 required.")
  } else if (minTSLength == 2 || minTSLength == 3) {
    return(2)
  } else {
    lowestEntropy <- Inf
    mOpt <- 2
    for (m in 2:min(7, minTSLength - 1)) {
      currentEntropy <- mean(sapply(tsList, pdcEntropy, m))
      if (currentEntropy < lowestEntropy) {
        mOpt <- m
        lowestEntropy <- currentEntropy
      }
    }
    return(mOpt)
  }
}

#' Multi-variate Permutation Distribution Entropy Heuristic
#'
#' Multi-variate version of \code{\link{pdcEntropyHeuristic}} which splits
#' multi-variate time series (matrices) in multiple uni-variate time series
#' (vectors) before applying the heuristic. This approach was proposed in
#' the dissertation of Brandmaier (2011, p. 25).
#'
#' @section References:
#'
#' Brandmaier, A. M. (2011). \emph{Permutation distribution clustering and structural
#' equation model trees} (Doctoral dissertation, Universität des Saarlandes, Saarbrücken).
#'
#' Brandmaier, A. M. (2015). pdc: An R Package for Complexity-Based Clustering of
#' Time Series. \emph{Journal of Statistical Software, 67}(5), 1-23.
#'
#' @param tsList A list of numeric vectors and/or matrixes which represent uni-variate
#' and/or multi-variate time series.
#' @return The optimal subsequence/pattern length as integer from the range [2,7].
#' @family PDC functions
#' @export
pdcEntropyHeuristicMult <- function(tsList) {
  tsList <- unlist(lapply(tsList, function(element) {
    if (is.matrix(element)) {
      lapply(seq_len(ncol(element)), function(j) element[,j])
    } else {
      return(list(element))
    }
  }), recursive = FALSE)
  return(pdcEntropyHeuristic(tsList))
}

#' Permutation Distribution Distance for Two Time Series
#'
#' Simplified version of the distance used for permutation distribution clustering
#' as described by Brandmaier (2015). Shortens the high-level code a bit for the two
#' time series case and handles some errors for short time series.
#'
#' @section References:
#'
#' Brandmaier, A. M. (2011). \emph{Permutation distribution clustering and structural
#' equation model trees} (Doctoral dissertation, Universität des Saarlandes, Saarbrücken).
#'
#' Brandmaier, A. M. (2015). pdc: An R Package for Complexity-Based Clustering of
#' Time Series. \emph{Journal of Statistical Software, 67}(5), 1-23.
#'
#' @importFrom pdc symmetricAlphaDivergence
#' @importFrom pdc codebook
#' @param x 1st numeric vector/time series (minimum length: 2 elements).
#' @param y 2nd numeric vector/time series (minimum length: 2 elements).
#' @param subSequenceLength Number of elements which form each subsequence. Will
#' be determined as integer in [2,7] by a heuristic if not provided.
#' @return The distance as double from the range [0,4].
#' @family PDC functions
#' @export
pdcDistTwoTS <- function(x, y, subSequenceLength = NULL) {
  if (is.null(subSequenceLength)) {
    subSequenceLength <- pdcEntropyHeuristic(list(x,y))
  }
  return(symmetricAlphaDivergence(codebook(x, m = subSequenceLength, t = 1),
                                  codebook(y, m = subSequenceLength, t = 1)))
}

#' Multi-variate Permutation Distribution Distance for Two Time Series
#'
#' Multi-variate version of \code{\link{pdcDistTwoTS}} which calculates the pdc
#' distance for each attribute and then takes the l2 norm of the resulting vector
#' as proposed by Brandmaier (2011, p. 18).
#'
#' @section References:
#'
#' Brandmaier, A. M. (2011). \emph{Permutation distribution clustering and structural
#' equation model trees} (Doctoral dissertation, Universität des Saarlandes, Saarbrücken).
#'
#' Brandmaier, A. M. (2015). pdc: An R Package for Complexity-Based Clustering of
#' Time Series. \emph{Journal of Statistical Software, 67}(5), 1-23.
#'
#' @importFrom pdc symmetricAlphaDivergence
#' @importFrom pdc codebook
#' @param x 1st numeric matrix (multi-variate time series; minimum number of rows: 2).
#' @param y 2nd numeric matrix (multi-variate time series; minimum number of rows: 2;
#' needs to have the same number of columns as the first time series).
#' @param subSequenceLength Number of elements which form each subsequence. Will
#' be determined as integer in [2,7] by a heuristic if not provided.
#' @return The distance as double from the range [0,4*attributeCount].
#' @family PDC functions
#' @export
pdcDistTwoTSMult <- function(x, y, subSequenceLength = NULL) {
  if (is.null(subSequenceLength)) {
    subSequenceLength <- pdcEntropyHeuristicMult(list(x,y))
  }
  # L2 norm of column distances
  return(sqrt(sum(sapply(seq_len(ncol(x)), function(j) {
    symmetricAlphaDivergence(codebook(x[, j], m = subSequenceLength, t = 1),
                              codebook(y[, j], m = subSequenceLength, t = 1))
  })^2)))
}

#' Permutation Distribution Distance for Two Lists of Time Series
#'
#' Version of \code{\link{pdcDistTSList}} which computes the PDC distance from
#' each time series of the first list to each time series of the second list.
#'
#' @importFrom pdc symmetricAlphaDivergence
#' @importFrom pdc codebook
#' @param tsList1 A list of numeric vectors (uni-variate time series). Each time
#' series should have at least two elements.
#' @param tsList2 A list of numeric vectors (uni-variate time series). Each time
#' series should have at least two elements. Is used for the entropy heuristic
#' if no subSequenceLength is provided,
#' @param subSequenceLength Number of elements which form each subsequence. Will
#' be determined as integer in [2,7] by a heuristic if not provided.
#' @return The distance matrix with each entry being a numeric from the range [0,4].
pdcDistTwoTSLists <- function(tsList1, tsList2, subSequenceLength = NULL) {
  if (is.null(subSequenceLength)) {
    # for comparability with previous distance computations which were conducted
    # with only one of the lists
    subSequenceLength <- pdcEntropyHeuristic(tsList2)
  }
  codebookList1 <- lapply(tsList1, codebook, m = subSequenceLength, t = 1)
  codebookList2 <- lapply(tsList2, codebook, m = subSequenceLength, t = 1)
  rowCount <- length(tsList1)
  colCount <- length(tsList2)
  distMat <- matrix(0, nrow = rowCount, ncol = colCount)
  for (i in 1:rowCount) {
    for (j in 1:colCount) {
      distMat[i,j] <- symmetricAlphaDivergence(codebookList1[[i]], codebookList2[[j]])
    }
  }
  return(distMat)
}

#' Permutation Distribution Distance for One List of Time Series
#'
#' Version of \code{\link{pdcDistTSList}} which computes the PDC distance from
#' each time series of a list to each time series of the same list.
#'
#' @importFrom pdc symmetricAlphaDivergence
#' @importFrom pdc codebook
#' @param tsList A list of numeric vectors (uni-variate time series). Each time
#' series should have at least two elements.
#' @param subSequenceLength Number of elements which form each subsequence. Will
#' be determined as integer in [2,7] by a heuristic if not provided.
#' @return The distance matrix with each entry being a numeric from the range [0,4].
pdcDistOneTSList <- function(tsList, subSequenceLength = NULL) {
  if (is.null(subSequenceLength)) {
    subSequenceLength <- pdcEntropyHeuristic(tsList)
  }
  codebookList <- lapply(tsList, codebook, m = subSequenceLength, t = 1)
  tsCount <- length(tsList)
  distMat <- matrix(0, nrow = tsCount, ncol = tsCount)
  for (i in 1:(tsCount - 1)) {
    for (j in (i + 1):tsCount) {
      distMat[i,j] <- symmetricAlphaDivergence(codebookList[[i]], codebookList[[j]])
      distMat[j,i] <- distMat[i,j]
    }
  }
  return(distMat)
}

#' Permutation Distribution Distance for a List of Time Series
#'
#' Simplified version of the distance used for permutation distribution clustering
#' as described by Brandmaier (2015). Shortens the high-level code a bit and
#' handles some errors for short time series.
#'
#' @section References:
#'
#' Brandmaier, A. M. (2011). \emph{Permutation distribution clustering and structural
#' equation model trees} (Doctoral dissertation, Universität des Saarlandes, Saarbrücken).
#'
#' Brandmaier, A. M. (2015). pdc: An R Package for Complexity-Based Clustering of
#' Time Series. \emph{Journal of Statistical Software, 67}(5), 1-23.
#'
#' @importFrom pdc symmetricAlphaDivergence
#' @importFrom pdc codebook
#' @param tsList 1) A list of numeric vectors (uni-variate time series). The
#' dissimilarity of the list to itself (each time series to each time series)
#' will be computed, resulting in a symmetric dissimilarity matrix. 2) A list
#' with two components, each being a list of numeric vectors (uni-variate time
#' series). The dissimilarity of each time series from the 1st component to each
#' time series from the 2nd component will be computed. The entropy heuristic is
#' only computed with the second component (for comparability to previously
#' computed distances). 1+2) Each time series should have at least two elements.
#' @param subSequenceLength Number of elements which form each subsequence. Will
#' be determined as integer in [2,7] by a heuristic if not provided.
#' @return The distance matrix with each entry being a numeric from the range [0,4].
#' @family PDC functions
#' @export
pdcDistTSList <- function(tsList, subSequenceLength = NULL) {
  if (length(tsList) == 2 && is.list(tsList[[1]]) && is.list(tsList[[2]])) {
    return(pdcDistTwoTSLists(tsList1 = tsList[[1]], tsList2 = tsList[[2]],
                             subSequenceLength = subSequenceLength))
  } else {
    return(pdcDistOneTSList(tsList = tsList,
                            subSequenceLength = subSequenceLength))
  }
}

#' Multi-variate Permutation Distribution Distance for Two Lists of Time Series
#'
#' Version of \code{\link{pdcDistTSListMult}} which computes the PDC distance from
#' each time series of the first list to each time series of the second list.
#'
#' @importFrom pdc symmetricAlphaDivergence
#' @importFrom pdc codebook
#' @param tsList1 A list of numeric matrixes (multi-variate time series). Each
#' time series should have at least two elements.
#' @param tsList2 A list of numeric metrixes (multi-variate time series). Each
#' time series should have at least two elements and the same number of attributes
#' as the series from the first list. Is used for the entropy heuristic if no
#' subSequenceLength is provided.
#' @param subSequenceLength Number of elements which form each subsequence. Will
#' be determined as integer in [2,7] by a heuristic if not provided.
#' @return The distance matrix with each entry being a numeric from the range
#' [0,4*attributeCount].
pdcDistTwoTSListsMult <- function(tsList1, tsList2, subSequenceLength = NULL) {
  rowCount <- length(tsList1)
  colCount <- length(tsList2)
  attributeCount <- unique(c(sapply(tsList1, ncol), sapply(tsList2, ncol)))
  if (length(attributeCount) > 1) {
    stop("All time series should have the same number of attributes.")
  }
  vectorTSList1 <- lapply(tsList1, function(element) {
    lapply(seq_len(attributeCount), function(j) element[,j])
  }) # list of time series transformed to list of list of single ts attributes
  vectorTSList2 <- lapply(tsList2, function(element) {
    lapply(seq_len(attributeCount), function(j) element[,j])
  })
  if (is.null(subSequenceLength)) {
    subSequenceLength <- pdcEntropyHeuristic(unlist(vectorTSList2, recursive = FALSE))
  }
  codebookList1 <- lapply(vectorTSList1, function(oneTSAttributes)
    lapply(oneTSAttributes, codebook, m = subSequenceLength, t = 1))
  codebookList2 <- lapply(vectorTSList2, function(oneTSAttributes)
    lapply(oneTSAttributes, codebook, m = subSequenceLength, t = 1))
  distMat <- matrix(0, nrow = rowCount, ncol = colCount)
  for (i in 1:rowCount) {
    for (j in 1:colCount) {
      distMat[i,j] <- sqrt(sum(sapply(seq_len(attributeCount), function(attributeNr) {
        symmetricAlphaDivergence(codebookList1[[i]][[attributeNr]],
                                 codebookList2[[j]][[attributeNr]])
      })^2))
    }
  }
  return(distMat)
}

#' Multi-variate Permutation Distribution Distance for One List of Time Series
#'
#' Version of \code{\link{pdcDistTSListMult}} which computes the PDC distance from
#' each time series of a list to each time series of the same list.
#'
#' @importFrom pdc symmetricAlphaDivergence
#' @importFrom pdc codebook
#' @param tsList A list of numeric matrixes (multi-variate time series). Each
#' time series should have at least two elements.
#' @param subSequenceLength Number of elements which form each subsequence. Will
#' be determined as integer in [2,7] by a heuristic if not provided.
#' @return The (symmetric) distance matrix with each entry being a numeric from
#' the range [0,4*attributeCount].
pdcDistOneTSListMult <- function(tsList, subSequenceLength = NULL) {
  tsCount <- length(tsList)
  attributeCount <- unique(sapply(tsList, ncol))
  if (length(attributeCount) > 1) {
    stop("All time series should have the same number of attributes.")
  }
  vectorTSList <- lapply(tsList, function(element) {
    lapply(seq_len(attributeCount), function(j) element[,j])
  }) # list of time series transformed to list of list of single ts attributes
  if (is.null(subSequenceLength)) {
    subSequenceLength <- pdcEntropyHeuristic(unlist(vectorTSList, recursive = FALSE))
  }
  codebookList <- lapply(vectorTSList, function(oneTSAttributes)
    lapply(oneTSAttributes, codebook, m = subSequenceLength, t = 1))
  distMat <- matrix(0, nrow = tsCount, ncol = tsCount)
  for (i in 1:(tsCount - 1)) {
    for (j in (i + 1):tsCount) {
      distMat[i,j] <- sqrt(sum(sapply(seq_len(attributeCount), function(attributeNr) {
        symmetricAlphaDivergence(codebookList[[i]][[attributeNr]],
                                 codebookList[[j]][[attributeNr]])
      })^2))
      distMat[j,i] <- distMat[i,j]
    }
  }
  return(distMat)
}

#' Multi-variate Permutation Distribution Distance for a List of Time Series
#'
#' Multi-variate version of \code{\link{pdcDistTSList}}.
#'
#' @section References:
#'
#' Brandmaier, A. M. (2011). \emph{Permutation distribution clustering and structural
#' equation model trees} (Doctoral dissertation, Universität des Saarlandes, Saarbrücken).
#'
#' Brandmaier, A. M. (2015). pdc: An R Package for Complexity-Based Clustering of
#' Time Series. \emph{Journal of Statistical Software, 67}(5), 1-23.
#'
#' @importFrom pdc symmetricAlphaDivergence
#' @importFrom pdc codebook
#' @param tsList 1) A list of numeric matrixes (multi-variate time series). The
#' dissimilarity of the list to itself (each time series to each time series)
#' will be computed, resulting in a symmetric dissimilarity matrix. 2) A list
#' with two components, each being a list of numeric matrixes (multi-variate time
#' series). The dissimilarity of each time series from the 1st component to each
#' time series from the 2nd component will be computed. The entropy heuristic is
#' only computed with the second component (for comparability to previously
#' computed distances). 1+2) Each time series should each have a minimum element
#' (row) count of two and all the same number of attributes (columns)).
#' @param subSequenceLength Number of elements which form each subsequence. Will
#' be determined as integer in [2,7] by a heuristic if not provided.
#' @return The distance matrix with each entry being a numeric from the range
#' [0,4*attributeCount].
#' @family PDC functions
#' @export
pdcDistTSListMult <- function(tsList, subSequenceLength = NULL) {
  if (length(tsList) == 2 && is.list(tsList[[1]]) && is.list(tsList[[2]])) {
    return(pdcDistTwoTSListsMult(tsList1 = tsList[[1]], tsList2 = tsList[[2]],
                                 subSequenceLength = subSequenceLength))
  } else {
    return(pdcDistOneTSListMult(tsList = tsList,
                                subSequenceLength = subSequenceLength))
  }
}

#' Compute cross-distance matrix
#'
#' Computes the cross-distance matrix between all (time series) objects from a
#' list, assuming a symmetric dissimilarity and a self-distance of 0. Uses
#' parallelization.
#'
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach  "%dopar%" foreach
#' @importFrom parallel makeCluster stopCluster
#' @param tsList A list containing the objects (e.g. numeric vectors representing
#' time series).
#' @param distMethod The name of a function which should be used to calculate the
#' distance between two objects.
#' @param distArgs Further parameters which should be passed to the distance
#' function in the form of a named list.
#' @param distExports Function (not part of a package) which are need for the
#' distance computation and therefore need to be exported to work with
#'  \code{foreach}.
#' @param distPackages Packages which are necessary for using the \code{distMethod}.
#' @param cpuCores The number of cores/threads which should be used for
#' parallelization. Setting this value to 1 forces sequential computation.
#' @param maxSequentialTime A short speed test is performed to decide if
#' parallelization is worth it (as it comes with a certain overhead). If the
#' speed test expects a total runtime <= the parameter value (in seconds), it
#' calls the distance method sequentially. A parameter value <= 0 forces parallel
#' execution (without speed test).
#' @param trace Print status information? (currently only the computation time
#' after finishing the computation).
#' @return A numeric matrix with row and column count == number of objects in
#' \code{tsList}, containing the pairwise distances.
#' @export
tsCrossDistMat <- function(tsList, distMethod, distArgs = list(), distExports = NULL,
                           distPackages = NULL, cpuCores = 2, maxSequentialTime = 3,
                           trace = TRUE) {
  startTime <- as.numeric(Sys.time())
  tsCount <- length(tsList)
  if (cpuCores > 1 && maxSequentialTime > 0) {# Speed test to determine if parallelization is worth it
    numSpeedTestSamples <- min(10, tsCount)
    speedTestSamples <- sample(1:tsCount, size = numSpeedTestSamples)
    speedTestStartTime <- as.numeric(Sys.time())
    for (i in 2:length(speedTestSamples)) {
      do.call(what = distMethod,
              args = c(list(tsList[[speedTestSamples[1]]], tsList[[speedTestSamples[i]]]), distArgs))
    }
    speedTestEndTime <- as.numeric(Sys.time())
    expectedComputationTime <- tsCount * (tsCount-1) / 2 / (numSpeedTestSamples-1) *
      (speedTestEndTime-speedTestStartTime)
  }
  # Real distance matrix computation
  if (cpuCores == 1 || (maxSequentialTime > 0 && expectedComputationTime <= maxSequentialTime)) {
    distMat <- sapply(1:tsCount, function(j) { # colum-wise
      if (j == 1) {
        return(c(0.0, rep(NA, times = tsCount-1)))
      } else {
        return(c(sapply(1:(j-1), function(i) {
          return(do.call(what = distMethod,
                         args = c(list(tsList[[i]], tsList[[j]]), distArgs)))
        }), 0.0, rep(NA, times = tsCount-j)))
      }
    })
  } else {
    computingCluster <- makeCluster(cpuCores)
    registerDoParallel(computingCluster)
    # Pairwise distance, but matrix assumed to be symmetric to main diagonal (which is 0)
    distMat <- rbind(foreach(i = 1:(tsCount-1), .export = distExports,
                             .packages = distPackages, .combine = "rbind") %dopar% {
      return(c(rep(NA, times = i-1), 0.0, sapply((i+1):tsCount, function(j) {
        return(do.call(what = distMethod, args = c(list(tsList[[i]], tsList[[j]]), distArgs)))
      })))
   }, c(rep(NA, times = tsCount-1), 0.0))
    stopCluster(computingCluster)
  }
  distMat[lower.tri(distMat, diag = FALSE)] <- t(distMat)[lower.tri(distMat, diag = FALSE)]
  endTime <- as.numeric(Sys.time())
  if (trace) {
    cat("Distance computation took ", (endTime - startTime), " secs.\n")
  }
  return(distMat)
}

#' Z-normalization
#'
#' Calculates the z-normalization by subtracting the mean and dividing by the
#' standard deviation (constant series are normalized to 0). Multi-variate time
#' series (matrices) are normalized per column.
#'
#' @importFrom stats sd
#' @param x A numeric vector/matrix (uni- or multivariate time series).
#' @return Z-normalized time series.
#' @export
znormalize <- function(x) {
  if (is.matrix(x)) {
    return(apply(x, 2, znormalize))
  }
  sdev <- sd(x)
  if (is.na(sdev) || sdev == 0) {
    return(rep(0, times = length(x)))
  } else {
    return((x - mean(x)) / sdev)
  }
}
