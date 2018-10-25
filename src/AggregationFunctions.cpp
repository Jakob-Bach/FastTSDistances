#include <Rcpp.h>
using namespace Rcpp;

//' Symbolic Aggregate Aproximation
//'
//' Converts a numeric time series to a series of discrete values (here: integers,
//' but could also be letters or other symbols), depending on where would they be
//' placed in the series of intervals called \code{limits}. Inspired by the pure R
//' implementation in \code{TSclust::convert.to.SAX.symbol()}, but faster because
//' C++ code.
//'
//' @section References:
//'
//' Lin, J., Keogh, E., Lonardi, S. & Chiu, B. (2003). A symbolic representation
//' of time series, with implications for streaming algorithms. In \emph{Proceedings
//' of the 8th acm sigmod workshop on research issues in data mining and knowledge
//' discovery} (pp. 2-11). ACM.
//'
//' @param x Numeric vector/time series.
//' @param limits The interval boundaries which determine the mapping to symbols.
//' A time series value in the i-th interval will be mapped to the value i. Your
//' limits should be a monotonically increasing vector, starting from negative
//' infinity and ending with positive infinity (length = number of symbols + 1).
//' @return The series converted to integer values.
//' @family SAX functions
//' @export
// [[Rcpp::export]]
NumericVector SAX_fast(NumericVector x, NumericVector limits) {
  NumericVector result = NumericVector(x.size());
  int symbol;
  for (int i = 0; i < x.size(); i++) {
    for(symbol = 1; x[i] > limits[symbol]; symbol++);
    result[i] = symbol;
  }
  return result;
}

//' Piecewise Aggregate Approximation
//'
//' Divides a time series into \code{windowCount} frames of equal length and
//' represents each interval by its mean. If the time series length is not
//' divisible by the number of windows, one element might belong to two windows,
//' but only with a certain fraction to each. Inspired by the pure R implementation
//' in \code{TSclust::PAA()}, but faster because C++ code.
//'
//' @section References:
//'
//' Keogh, E. J. & Pazzani, M. J. (2000). Scaling up dynamic time warping for
//' datamining applications. In \emph{Proceedings of the sixth acm sigkdd
//' international conference on knowledge discovery and data mining} (pp. 285-289).
//' ACM.
//'
//' Keogh, E., Chakrabarti, K., Pazzani, M. & Mehrotra, S. (2001). Dimensionality
//' reduction for fast similarity search in large time series databases.
//' \emph{Knowledge and information Systems}, 3(3), 263-286.
//'
//' @param x Numeric vector/time series.
//' @param windowCount The number of windows for the shortened time series.
//' @return A numeric vector of the length \code{windowCount}.
//' @family piecewise aggregation functions
//' @export
// [[Rcpp::export]]
NumericVector PAA_fast(NumericVector x, int windowCount) {
  NumericVector result = NumericVector(windowCount);
  double windowLength = x.size() * 1.0 / windowCount;
  // Indices which do not fall fully into window (e.g. 5 and 8 if window goes from
  // 5.7 to 7.4) are considered proportionally
  double windowStartIdx, windowEndIdx, startFraction, endFraction, aggValue;
  int fullWindowStartIdx, fullWindowEndIdx;
  for (int i = 0; i < windowCount; i++) {
    windowStartIdx = windowLength * i;
    windowEndIdx = windowLength * (i+1) - 1;
    fullWindowStartIdx = std::ceil(windowStartIdx);
    fullWindowEndIdx = std::floor(windowEndIdx);
    startFraction = fullWindowStartIdx - windowStartIdx;
    endFraction = windowEndIdx - fullWindowEndIdx;
    aggValue = 0;
    if (startFraction > 0) {
      aggValue += x[fullWindowStartIdx - 1] * startFraction;
    }
    if (endFraction > 0) {
      aggValue += x[fullWindowEndIdx + 1] * endFraction;
    }
    for (int j = fullWindowStartIdx; j <= fullWindowEndIdx; j++) {
      aggValue += x[j];
    }
    result[i] = aggValue / windowLength;
  }
  return result;
}

//' Piecewise Maximum Aggregate Approximation
//'
//' Divides a time series into \code{windowCount} frames of equal length and
//' represents each interval by its maximum. If the time series length is not
//' divisible by the number of windows, elements still assigned uniquely to one
//' window each.
//'
//' The function is similar to \code{\link{PAA_fast}}, but simply uses a different
//' aggregate.
//'
//' @param x Numeric vector/time series.
//' @param windowCount The number of windows for the shortened time series.
//' @return A numeric vector of the length \code{windowCount}.
//' @family piecewise aggregation functions
//' @export
// [[Rcpp::export]]
NumericVector PMaxAA_fast(NumericVector x, int windowCount) {
  NumericVector result = NumericVector(windowCount);
  double windowLength = x.size() * 1.0 / windowCount;
  double aggValue;
  int windowStartIdx, windowEndIdx;
  for (int i = 0; i < windowCount; i++) {
    windowStartIdx = std::floor(windowLength * i);
    windowEndIdx = std::floor(windowLength * (i+1) - 1);
    aggValue = x[windowStartIdx];
    for (int j = windowStartIdx + 1; j <= windowEndIdx; j++) {
      aggValue = std::max(aggValue, x[j]);
    }
    result[i] = aggValue;
  }
  return result;
}

//' Piecewise Median Aggregate Approximation
//'
//' Divides a time series into \code{windowCount} frames of equal length and
//' represents each interval by its median. If the time series length is not
//' divisible by the number of windows, elements still assigned uniquely to one
//' window each.
//'
//' The function is similar to \code{\link{PAA_fast}}, but simply uses a different
//' aggregate.
//'
//' @param x Numeric vector/time series.
//' @param windowCount The number of windows for the shortened time series.
//' @return A numeric vector of the length \code{windowCount}.
//' @family piecewise aggregation functions
//' @export
// [[Rcpp::export]]
NumericVector PMedAA_fast(NumericVector x, int windowCount) {
  NumericVector result = NumericVector(windowCount);
  double windowLength = x.size() * 1.0 / windowCount;
  // Indices which do not fall fully into window (e.g. 5 and 8 if window goes from
  // 5.7 to 7.4) are considered proportionally
  double med1, med2, medIndexDouble;
  int windowStartIdx, windowEndIdx, medIndexInt;
  for (int i = 0; i < windowCount; i++) {
    windowStartIdx = std::floor(windowLength * i);
    windowEndIdx = std::floor(windowLength * (i+1) - 1);
    // median code adapted from http://gallery.rcpp.org/articles/robust-estimators/
    NumericVector window = x[Range(windowStartIdx, windowEndIdx)];
    medIndexDouble = (window.size() - 1) / 2.0;
    medIndexInt = (int) medIndexDouble;
    if (medIndexDouble == medIndexInt) {
      std::nth_element(window.begin(), window.begin() + medIndexInt, window.end());
      result[i] = window[medIndexInt];
    } else {
      std::nth_element(window.begin(), window.begin() + medIndexInt + 1, window.end());
      med1 = window[medIndexInt + 1];
      std::nth_element(window.begin(), window.begin() + medIndexInt,
                       window.begin() + medIndexInt + 1);
      med2 = window[medIndexInt];
      result[i] = 0.5 * (med1 + med2);
    }
  }
  return result;
}

//' Piecewise Minimum Aggregate Approximation
//'
//' Divides a time series into \code{windowCount} frames of equal length and
//' represents each interval by its minimum. If the time series length is not
//' divisible by the number of windows, elements still assigned uniquely to one
//' window each.
//'
//' The function is similar to \code{\link{PAA_fast}}, but simply uses a different
//' aggregate.
//'
//' @param x Numeric vector/time series.
//' @param windowCount The number of windows for the shortened time series.
//' @return A numeric vector of the length \code{windowCount}.
//' @family piecewise aggregation functions
//' @export
// [[Rcpp::export]]
NumericVector PMinAA_fast(NumericVector x, int windowCount) {
  NumericVector result = NumericVector(windowCount);
  double windowLength = x.size() * 1.0 / windowCount;
  double aggValue;
  int windowStartIdx, windowEndIdx;
  for (int i = 0; i < windowCount; i++) {
    windowStartIdx = std::floor(windowLength * i);
    windowEndIdx = std::floor(windowLength * (i+1) - 1);
    aggValue = x[windowStartIdx];
    for (int j = windowStartIdx + 1; j <= windowEndIdx; j++) {
      aggValue = std::min(aggValue, x[j]);
    }
    result[i] = aggValue;
  }
  return result;
}

//' Piecewise Standard Deviation Aggregate Approximation
//'
//' Divides a time series into \code{windowCount} frames of equal length and
//' represents each interval by its standard deviation. If the time series length
//' is not divisible by the number of windows, one element might belong to two
//' windows, but only with a certain fraction to each.
//'
//' The function is similar to \code{\link{PAA_fast}}, but simply uses a different
//' aggregate.
//'
//' @param x Numeric vector/time series.
//' @param windowCount The number of windows for the shortened time series.
//' @param sample Compute sample standard deviation instead of population
//' standard deviation (divide by n-1 instead of n).
//' @return A numeric vector of the length \code{windowCount}.
//' @family piecewise aggregation functions
//' @export
// [[Rcpp::export]]
NumericVector PSDAA_fast(NumericVector x, int windowCount, bool sample = false) {
  NumericVector result = NumericVector(windowCount);
  double windowLength = x.size() * 1.0 / windowCount;
  // Indices which do not fall fully into window (e.g. 5 and 8 if window goes from
  // 5.7 to 7.4) are considered proportionally
  double windowStartIdx, windowEndIdx, startFraction, endFraction, linearSum, squaredSum;
  int fullWindowStartIdx, fullWindowEndIdx;
  for (int i = 0; i < windowCount; i++) {
    windowStartIdx = windowLength * i;
    windowEndIdx = windowLength * (i+1) - 1;
    fullWindowStartIdx = std::ceil(windowStartIdx);
    fullWindowEndIdx = std::floor(windowEndIdx);
    startFraction = fullWindowStartIdx - windowStartIdx;
    endFraction = windowEndIdx - fullWindowEndIdx;
    linearSum = 0;
    squaredSum = 0;
    if (startFraction > 0) {
      linearSum = x[fullWindowStartIdx - 1] * startFraction;
      squaredSum = x[fullWindowStartIdx - 1] * x[fullWindowStartIdx - 1] * startFraction;
    }
    if (endFraction > 0) {
      linearSum += x[fullWindowEndIdx + 1] * endFraction;
      squaredSum += x[fullWindowEndIdx + 1] * x[fullWindowEndIdx + 1] * endFraction;
    }
    for (int j = fullWindowStartIdx; j <= fullWindowEndIdx; j++) {
      linearSum += x[j];
      squaredSum += x[j] * x[j];
    }
    result[i] = fmax((squaredSum / windowLength - pow(linearSum / windowLength, 2)), 0);
    if (sample) {
      result[i] = std::sqrt(windowLength * result[i] / (windowLength - 1));
    } else {
      result[i] = std::sqrt(result[i]);
    }
  }
  return result;
}

//' Mean with Weighted First And Last Element
//'
//' Calculates the mean of a sub-vector, weighting start and end if the corresponding
//' indices are no whole numbers.
//'
//' @param x Numeric vector.
//' @param startIdx Index of first element of sub-vector.
//' @param endIdx Index of last element of sub-vector.
//' @return The weighted mean as double.
// [[Rcpp::export]]
double subVectorMean_fast(NumericVector x, double startIdx, double endIdx) {
  double fullWindowStartIdx = std::ceil(startIdx);
  double fullWindowEndIdx = std::floor(endIdx);
  double startFraction = fullWindowStartIdx - startIdx;
  double endFraction = fullWindowEndIdx - endIdx;
  double result = 0;
  if (startFraction > 0) {
    result += x[fullWindowStartIdx - 1] * startFraction;
  }
  if (endFraction > 0) {
    result += x[fullWindowEndIdx + 1] * endFraction;
  }
  for (int j = fullWindowStartIdx; j <= fullWindowEndIdx; j++) {
    result += x[j];
  }
  return result / (endIdx - startIdx + 1);
}

//' Piecewise Skewness Aggregate Approximation
//'
//' Divides a time series into \code{windowCount} frames of equal length and
//' represents each interval by its skewness. If the time series length is not
//' divisible by the number of windows, one element might belong to two windows,
//' but only with a certain fraction to each.
//'
//' The function is similar to \code{\link{PAA_fast}}, but simply uses a different
//' aggregate.
//'
//' @param x Numeric vector/time series.
//' @param windowCount The number of windows for the shortened time series.
//' @param nanReplace All NaN elements will be replaced with this value (default: NaN).
//' NaNs can occur if a segment is constant.
//' @return A numeric vector of the length \code{windowCount}.
//' @family piecewise aggregation functions
//' @export
// [[Rcpp::export]]
NumericVector PSkewAA_fast(NumericVector x, int windowCount, double nanReplace = NA_REAL) {
  if (std::isnan(nanReplace)) { // replace with proper NaN
    nanReplace = std::numeric_limits<double>::quiet_NaN();
  }
  NumericVector result = NumericVector(windowCount);
  double windowLength = x.size() * 1.0 / windowCount;
  // Indices which do not fall fully into window (e.g. 5 and 8 if window goes from
  // 5.7 to 7.4) are considered proportionally
  double windowStartIdx, windowEndIdx, startFraction, endFraction, windowMean,
  windowSD, squaredSum, cubedSum;
  int fullWindowStartIdx, fullWindowEndIdx;
  for (int i = 0; i < windowCount; i++) {
    windowStartIdx = windowLength * i;
    windowEndIdx = windowLength * (i+1) - 1;
    fullWindowStartIdx = std::ceil(windowStartIdx);
    fullWindowEndIdx = std::floor(windowEndIdx);
    startFraction = fullWindowStartIdx - windowStartIdx;
    endFraction = windowEndIdx - fullWindowEndIdx;
    windowMean = subVectorMean_fast(x, windowStartIdx, windowEndIdx);
    squaredSum = 0;
    cubedSum = 0;
    if (startFraction > 0) {
      squaredSum += pow(x[fullWindowStartIdx - 1], 2) * startFraction;
      cubedSum += pow(x[fullWindowStartIdx - 1] - windowMean, 3) * startFraction;
    }
    if (endFraction > 0) {
      squaredSum += pow(x[fullWindowEndIdx + 1], 2) * endFraction;
      cubedSum += pow(x[fullWindowEndIdx + 1] - windowMean, 3) * endFraction;
    }
    for (int j = fullWindowStartIdx; j <= fullWindowEndIdx; j++) {
      squaredSum += x[j] * x[j];
      cubedSum += pow(x[j] - windowMean, 3);
    }
    windowSD = sqrt(fmax((squaredSum / windowLength - pow(windowMean, 2)), 0));
    result[i] = cubedSum / (windowLength * pow(windowSD, 3));
    // Alternative checks for constant series (depending on numerics, architecture)
    if (windowSD == 0 || std::isnan(result[i]) || std::isinf(result[i])) {
      result[i] = nanReplace;
    }
  }
  return result;
}

//' Piecewise Kurtosis Aggregate Approximation
//'
//' Divides a time series into \code{windowCount} frames of equal length and
//' represents each interval by its kurtosis. If the time series length is not
//' divisible by the number of windows, one element might belong to two windows,
//' but only with a certain fraction to each.
//'
//' The function is similar to \code{\link{PAA_fast}}, but simply uses a different
//' aggregate.
//'
//' @param x Numeric vector/time series.
//' @param windowCount The number of windows for the shortened time series.
//' @param nanReplace All NaN elements will be replaced with this value (default: NaN).
//' NaNs can occur if a segment is constant.
//' @param excess Compute excess kurtosis? (kurtosis - 3, is zero for normal distribution)
//' @return A numeric vector of the length \code{windowCount}.
//' @family piecewise aggregation functions
//' @export
// [[Rcpp::export]]
NumericVector PKurtAA_fast(NumericVector x, int windowCount, double nanReplace = NA_REAL,
                           bool excess = false) {
  if (std::isnan(nanReplace)) { // replace with proper NaN
    nanReplace = std::numeric_limits<double>::quiet_NaN();
  }
  NumericVector result = NumericVector(windowCount);
  double windowLength = x.size() * 1.0 / windowCount;
  // Indices which do not fall fully into window (e.g. 5 and 8 if window goes from
  // 5.7 to 7.4) are considered proportionally
  double windowStartIdx, windowEndIdx, startFraction, endFraction, windowMean,
  windowVar, squaredSum, quadSum;
  int fullWindowStartIdx, fullWindowEndIdx;
  for (int i = 0; i < windowCount; i++) {
    windowStartIdx = windowLength * i;
    windowEndIdx = windowLength * (i+1) - 1;
    fullWindowStartIdx = std::ceil(windowStartIdx);
    fullWindowEndIdx = std::floor(windowEndIdx);
    startFraction = fullWindowStartIdx - windowStartIdx;
    endFraction = windowEndIdx - fullWindowEndIdx;
    windowMean = subVectorMean_fast(x, windowStartIdx, windowEndIdx);
    squaredSum = 0;
    quadSum = 0;
    if (startFraction > 0) {
      squaredSum += pow(x[fullWindowStartIdx - 1], 2) * startFraction;
      quadSum += pow(x[fullWindowStartIdx - 1] - windowMean, 4) * startFraction;
    }
    if (endFraction > 0) {
      squaredSum += pow(x[fullWindowEndIdx + 1], 2) * endFraction;
      quadSum += pow(x[fullWindowEndIdx + 1] - windowMean, 4) * endFraction;
    }
    for (int j = fullWindowStartIdx; j <= fullWindowEndIdx; j++) {
      squaredSum += x[j] * x[j];
      quadSum += pow(x[j] - windowMean, 4);
    }
    windowVar = fmax(squaredSum / windowLength - pow(windowMean, 2), 0);
    result[i] = quadSum / (windowLength * windowVar * windowVar);
    // Alternative checks for constant series (depending on numerics, architecture)
    if (windowVar == 0 || std::isnan(result[i]) || std::isinf(result[i])) {
      result[i] = nanReplace;
    }
  }
  if (excess) {
    result = result - 3;
  }
  return result;
}
