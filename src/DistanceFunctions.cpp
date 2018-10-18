#include <Rcpp.h>
using namespace Rcpp;

//' Average Univariate Time Series
//'
//' Multiplies a list of univariate time series with a weight vector and
//' sums up the result to get a weighted average.
//'
//' @param tsList A list of NumericVectors all having the same length.
//' @param weights A vector of weights, having the same length as the list.
//' @return The weighted average of the series from the list.
//' @export
// [[Rcpp::export]]
NumericVector averageTimeSeries_fast(List tsList, NumericVector weights) {
  NumericVector result = NumericVector(((NumericVector) tsList[0]).size());
  for (int i = 0; i < weights.size(); i++) {
    result = result + ((NumericVector) tsList[i]) * weights[i];
  }
  return result;
}

//' Average Multivariate Time Series
//'
//' Multiplies a list of multivariate time series with a weight vector and
//' sums up the result to get a weighted average.
//'
//' @param tsList A list of NumericMatrix all having the same length.
//' @param weights A vector of weights, having the same length as the list.
//' @return The weighted average of the series from the list.
//' @export
// [[Rcpp::export]]
NumericMatrix averageTimeSeriesMult_fast(List tsList, NumericVector weights) {
  NumericMatrix curMatrix = tsList[0];
  NumericMatrix result = NumericMatrix(curMatrix.rows(), curMatrix.cols());
  for (int i = 0; i < weights.size(); i++) {
    NumericMatrix curMatrix = tsList[i];
    for (int j = 0; j < result.ncol(); j++) {
      result(_,j) = result.column(j) + curMatrix.column(j) * weights[i];
    }
  }
  return result;
}

//' L2 Complexity Correction Factor for a Time Series Distance
//'
//' Calculates the complexity correction factor for the distance between two time
//' series, using the L2 norm of each time series' diff vector as complexity
//' estimate. Can be combined with any distance as a scaling factor (distances
//' between vectors of different complexity become more prominent). Does not obey
//' the triangular equality if its combined with the Euclidean distance, but a
//' relaxed version (see reference).
//'
//' This factor is currently integrated as a parameter into the L2 distance and
//' dynamic time warping distance of this package.
//'
//' @section References:
//' Batista, G. E., Keogh, E. J., Tataw, O. M. & De Souza, V. M. (2014). Cid: An
//' efficient complexity-invariant distance for time series. \emph{Data Mining and
//' Knowledge Discovery, 28}(3), 634-669.
//'
//' @param x 1st numeric vector/time series.
//' @param y 2nd numeric vector/time series.
//' @return The complexity correction factor as double. Is infinity if one series
//' is constant and the other one not.
//' @seealso \code{\link{l2Dist_fast}}, \code{\link{DTWDist_fast}}
//' @export
// [[Rcpp::export]]
double l2CompCorFactor_fast(NumericVector x, NumericVector y) {
  double xComplexity = sqrt(sum(diff(x)*diff(x)));
  double yComplexity = sqrt(sum(diff(y)*diff(y)));
  if (xComplexity == 0 && yComplexity == 0) {
    return 1; // both constant, so both same complexity, prevent 0/0
  } else {
    return fmax(xComplexity,yComplexity) / fmin(xComplexity, yComplexity);
  }
}

//' Multi-variate L2 Complexity Correction Factor
//'
//' Calculates the complexity correction factor for the distance between two
//' multi-variate time series, using the L2 norm of each time series' attributes'
//' diff vectors as complexity estimate. Can be combined with any distance as a
//' scaling factor (distances between vectors of different complexity become more
//' prominent). Does not obey the triangular equality if its combined with the
//' Euclidean distance, but a relaxed version (see reference).
//'
//' This factor is currently integrated as a parameter into the L2 distance and
//' dynamic time warping distance of this package.
//'
//' @section References:
//'
//' Batista, G. E., Keogh, E. J., Tataw, O. M. & De Souza, V. M. (2014). Cid: An
//' efficient complexity-invariant distance for time series. \emph{Data Mining and
//' Knowledge Discovery, 28}(3), 634-669.
//'
//' Kotsakos, D., Trajcevski, G., Gunopulos, D. & Aggarwal, C. C. (2014). Time-series
//' data clustering. In C. C. Aggarwal & C. K. Reddy (Eds.), \emph{Data clustering :
//' Algorithms and applications} (pp. 357–380). Chapman & Hall/CRC data mining and
//' knowledge discovery series. Boca Raton: CRC Press.
//'
//' @param x 1st numeric matrix/multi-variate time series.
//' @param y 2nd numeric vector/multi-variate time series.
//' @return The complexity correction factor as double. Is infinity if one series
//' is constant in all attributes and the other one not.
//' @seealso \code{\link{l2Dist_fast}}, \code{\link{DTWDist_fast}}
//' @export
// [[Rcpp::export]]
double l2CompCorFactorMult_fast(NumericMatrix x, NumericMatrix y) {
  double xComplexity = 0, yComplexity = 0;
  for (int j = 0; j < x.ncol(); j++) {
    xComplexity += sum(diff(x.column(j)) * diff(x.column(j)));
  }
  xComplexity = sqrt(xComplexity);
  for (int j = 0; j < y.ncol(); j++) {
    yComplexity += sum(diff(y.column(j)) * diff(y.column(j)));
  }
  yComplexity = sqrt(yComplexity);
  if (xComplexity == 0 && yComplexity == 0) {
    return 1; // both constant, so both same complexity, prevent 0/0
  } else {
    return fmax(xComplexity, yComplexity) / fmin(xComplexity, yComplexity);
  }
}

//' (Fast) L2 Norm
//'
//' Computes the standard Euclidean norm with a fast C++ implementation.
//'
//' @param x A numeric vector/time series.
//' @return The norm as double.
//' @family L_p distances
//' @export
// [[Rcpp::export]]
double l2Norm_fast(NumericVector x) {
  double result = 0;
  int n = x.size();
  for (int i = 0; i < n; i++) {
    result += x[i]*x[i];
  }
  return sqrt(result);
}

//' Temporal Correlation
//'
//' Calculates the temporal correlation (correlation between the difference vectors)
//' as defined by Chouakria and Nagabhushan (2007). We additionally set the
//' correlation between two constant series to 1 and between a constant and a
//' non-constant one to 0 to guarantee that the result is always a proper (aka
//' not NaN, not infinity) number.
//'
//' @section References:
//'
//' Chouakria, A. D. & Nagabhushan, P. N. (2007). Adaptive dissimilarity index
//' for measuring time series proximity. \emph{Advances in Data Analysis and
//' Classification, 1}(1), 5-21.
//'
//' @param x 1st numeric vector/time series.
//' @param y 2nd numeric vector/time series.
//' @return The temporal correlation as double from the range [-1,1].
// [[Rcpp::export]]
double tempCor_fast(NumericVector x, NumericVector y) {
  double lengthDiffX = l2Norm_fast(diff(x));
  double lengthDiffY = l2Norm_fast(diff(y));
  if (lengthDiffX == 0) { // handle division by zero cases
    if (lengthDiffY == 0) { // if both constant
      return 1;
    }
  } else if (lengthDiffY != 0) { // ordinary case
    return sum(diff(x)*diff(y)) / (lengthDiffX * lengthDiffY);
  }
  return 0; // if one constant, other one not
}

//' Temporal Correlation-Based Correction Factor for a Time Series Distance
//'
//' Considers the dissimilarity of two time series regarding their behavior,
//' namely if they move in the same direction (diff vectors used) at the
//' different points in time. This correlation between [-1,1] is scaled with
//' the exponential function to (0,2) (depending on \code{k}) and should be
//' multiplied with another dissimilarity to enhance it with this behavioral
//' information. Introduced by  Chouakria and Nagabhushan (2007).
//'
//' This factor is currently integrated as a parameter into the L2 distance and
//' dynamic time warping distance of this package.
//'
//' @section References:
//'
//' Chouakria, A. D. & Nagabhushan, P. N. (2007). Adaptive dissimilarity index
//' for measuring time series proximity. \emph{Advances in Data Analysis and
//' Classification, 1}(1), 5-21.
//'
//' @param x 1st numeric vector/time series.
//' @param y 2nd numeric vector/time series.
//' @param k A non-negative constant for scaling.
//' @return The correlation-based scaling factor as double from the range (0,2).
//' @export
// [[Rcpp::export]]
double cortFactor_fast(NumericVector x, NumericVector y, double k = 2) {
  return 2 / (1 + exp(k * tempCor_fast(x, y)));
}

//' Multi-variate Temporal Correlation-Based Correction Factor
//'
//' Considers the dissimilarity of two time series regarding their behavior,
//' namely if they move in the same direction (diff vectors used) at the
//' different points in time. This correlation between [-1,1] is firstly
//' calculated for each dimension/attribute separately and then averaged.
//' It is scaled with the exponential function to (0,2) (depending on \code{k})
//' and should be multiplied with another dissimilarity to enhance it with
//' this behavioral information. Introduced by Chouakria and Nagabhushan (2007)
//' for the uni-variate case.
//'
//' This factor is currently integrated as a parameter into the L2 distance and
//' dynamic time warping distance of this package.
//'
//' @section References:
//'
//' Chouakria, A. D. & Nagabhushan, P. N. (2007). Adaptive dissimilarity index
//' for measuring time series proximity. \emph{Advances in Data Analysis and
//' Classification, 1}(1), 5-21.
//'
//' Kotsakos, D., Trajcevski, G., Gunopulos, D. & Aggarwal, C. C. (2014). Time-series
//' data clustering. In C. C. Aggarwal & C. K. Reddy (Eds.), \emph{Data clustering :
//' Algorithms and applications} (pp. 357–380). Chapman & Hall/CRC data mining and
//' knowledge discovery series. Boca Raton: CRC Press.
//'
//' @param x 1st numeric matrix/multi-variate time series.
//' @param y 2nd numeric matrix/multi-variate time series.
//' @param k A non-negative constant for scaling.
//' @return The correlation-based scaling factor as double from the range (0,2).
//' @export
// [[Rcpp::export]]
double cortFactorMult_fast(NumericMatrix x, NumericMatrix y, double k = 2) {
  double avgTempCorrelation = 0;
  for (int j = 0; j < x.ncol(); j++) {
    avgTempCorrelation += tempCor_fast(x.column(j), y.column(j));
  }
  avgTempCorrelation /= x.ncol();
  return 2 / (1 + exp(k * avgTempCorrelation));
}

//' Pairwise Absolute Distance
//'
//' Computes the pairwise absolute distance between two numeric vectors
//' (e.g. univariate time series)
//'
//' @param x 1st numeric vector/time series.
//' @param y 2nd numeric vector/time series.
//' @return An |x|*|y| matrix containing the pairwise absolute distances.
//' @export
// [[Rcpp::export]]
NumericMatrix vectorCrossDistMat(NumericVector x, NumericVector y) {
  NumericMatrix distMat(x.size(), y.size());
  for (int i = 0; i < x.size(); i++) {
    for (int j = 0; j < y.size(); j++) {
      distMat(i,j) = fabs(x[i] - y[j]);
    }
  }
  return distMat;
}

//' (Fast) Correlation-based Dissimilarity
//'
//' Computes correlation-based dissimilarity as described by Golay et al. (1998).
//' The coding is inspired by the \code{TSclust::diss.cor()} method, but faster
//' because of the the C++ implementation.
//'
//' @section References:
//'
//' Golay, X., Kollias, S., Stoll, G., Meier, D., Valavanis, A. & Boesiger, P.
//' (1998). A new correlation-based fuzzy logic clustering algorithm for fmri.
//' \emph{Magnetic Resonance in Medicine, 40}(2), 249-260.
//'
//' @param x 1st numeric vector/time series.
//' @param y 2nd numeric vector/time series.
//' @param beta If this parameter is smaller/equal zero, the formula \eqn{dist(x,y)=
//' \sqrt{2*(1-cor(x,y))}}{dist(x,y) = sqrt(2*(1-cor(x,y)))} is used (d1 in the
//' paper, equals \eqn{\frac{l2Dist(x,y)}{\sqrt{n}}}{l2Dist(x,y)/sqrt(n)} if time
//' series are z-standardized), otherwise
//' \eqn{dist(x,y)=\sqrt{(\frac{1-cor(x,y)}{1+cor(x,y)})^{\beta}}}{dist(x,y) =
//' sqrt(((1-cor(x,y))/(1+cor(x,y)))^beta)} (called d2 in the paper).
//' @return The dissimilarity as double in the range [0,sqrt(2)] if \code{beta == 0}
//' and [0,Inf] otherwise. Is NaN if at least one series is constant.
//' @export
// [[Rcpp::export]]
double corDist_fast(NumericVector x, NumericVector y, double beta = 0) {
  double sum_x = 0, sum_x2 = 0, sum_y = 0, sum_y2 = 0, sum_xy = 0;
  int n = x.size();
  for (int i = 0; i < n; i++) {
    sum_x += x[i];
    sum_x2 += x[i]*x[i];
    sum_y += y[i];
    sum_y2 += y[i]*y[i];
    sum_xy += x[i]*y[i];
  }
  double cor = (n*sum_xy - sum_x*sum_y) /
    sqrt((n*sum_x2 - sum_x*sum_x) * (n*sum_y2 - sum_y*sum_y));
  if (std::isnan(cor)) { // one series constant
    return cor;
  } else {
    cor = fmin(fmax(cor, -1), 1); // else we could run into numerical problems
    if (beta <= 0) { // d1, no beta
      return sqrt(2*(1-cor));
    } else { // d2, beta used
      return sqrt(pow((1-cor)/(1+cor), beta));
    }
  }
}

//' (Fast) Edit Distance on Real Sequence
//'
//' Computes the Edit distance on Real Sequence as described by Chen, Özsu
//' and Oria (2005). A match between two (real-valued) time series elements exists
//' if their L1 distance is below an \code{epsilon}. Apart from that, the
//' computation is similar to the standard edit distance. The coding is inspired
//' by the \code{TSdist::EDRDistance()} method, but faster because point-to-point
//' distances computation is integrated into the C++ code.
//'
//' Despite the name, it is not really a distance in the strict sense, as EDR
//' violates the triangular inequality.
//'
//' @section References:
//' Chen, L., Özsu, M. T. & Oria, V. (2005). Robust and fast similarity search
//' for moving object trajectories. In \emph{Proceedings of the 2005 acm sigmod
//' international conference on management of data} (pp. 491–502). ACM.
//'
//' @param x 1st numeric vector/time series.
//' @param y 2nd numeric vector/time series.
//' @param epsilon Maximum distance between two time series elements to count a
//' match.
//' @param normalize Normalize the result to [0,1] considering the maximum
//' possible dissimilarity.
//' @return The distance as double.
//' @family Edit distance functions
//' @export
// [[Rcpp::export]]
double EDRDist_fast(NumericVector x, NumericVector y, double epsilon,
                    bool normalize = false) {
  // Create cost matrix; is computed row by row and we only need access to
  // the current plus the previous row, so two rows are sufficient
  IntegerMatrix costMat(2, y.size() + 1);
  // Fill 1st row
  costMat(0,0) = 0;
  for (int j = 1; j <= y.size(); j++) {
    costMat(0,j) = j;
  }
  // Find optimal path
  int bottom, top; // we don't shift values when changing rows, just change indices
  for (int i = 1; i <= x.size(); i++) {
    bottom = i % 2; // current row
    top = (i - 1) % 2; // previous row
    costMat(bottom, 0) = i;
    for (int j = 1; j <= y.size(); j++) {
      costMat(bottom, j) = fmin(fmin(costMat(top, j) + 1, costMat(bottom, j-1) + 1),
                                costMat(top, j-1) + (fabs(x[i-1] - y[j-1]) > epsilon));
    }
  }
  if (normalize) {
    return costMat(x.size() % 2, y.size()) / fmax(x.size(), y.size());
  } else {
    return costMat(x.size() % 2, y.size());
  }
}

//' (Fast) Multi-variate Edit Distance on Real Sequence
//'
//' Computes the Edit distance on Real Sequence as described by Chen, Özsu
//' and Oria (2005). A match between two time series elements exists
//' if the L1 distances between corresponding attributes are all is below an
//' \code{epsilon}. Apart from that, the computation is similar to the standard
//' edit distance. The coding is inspired by the \code{TSdist::EDRDistance()}
//' method, but faster because point-to-point distances computation is integrated
//' into the C++ code.
//'
//' Despite the name, it is not really a distance in the strict sense, as EDR
//' violates the triangular inequality.
//'
//' @section References:
//' Chen, L., Özsu, M. T. & Oria, V. (2005). Robust and fast similarity search
//' for moving object trajectories. In \emph{Proceedings of the 2005 acm sigmod
//' international conference on management of data} (pp. 491–502). ACM.
//'
//' @param x 1st numeric matrix/multi-variate time series.
//' @param y 2nd numeric matrix/multi-variate time series.
//' @param epsilon Maximum distance between two time series elements to count a
//' match.
//' @param normalize Normalize the result to [0,1] considering the maximum
//' possible dissimilarity.
//' @return The distance as double.
//' @family Edit distance functions
//' @export
// [[Rcpp::export]]
double EDRDistMult_fast(NumericMatrix x, NumericMatrix y, double epsilon,
                        bool normalize = false) {
  // Create cost matrix; is computed row by row and we only need access to
  // the current plus the previous row, so two rows are sufficient
  IntegerMatrix costMat(2, y.nrow() + 1);
  // Fill 1st row
  costMat(0,0) = 0;
  for (int j = 1; j <= y.nrow(); j++) {
    costMat(0,j) = j;
  }
  // Find optimal path
  int bottom, top; // we don't shift values when changing rows, just change indices
  for (int i = 1; i <= x.nrow(); i++) {
    bottom = i % 2; // current row in cost matrix
    top = (i - 1) % 2; // previous row in cost matrix
    costMat(bottom, 0) = i;
    for (int j = 1; j <= y.nrow(); j++) {
      costMat(bottom, j) = fmin(fmin(costMat(top, j) + 1, costMat(bottom, j-1) + 1),
              costMat(top, j-1) + any(abs(x.row(i-1) - y.row(j-1)) > epsilon).is_true());
    }
  }
  if (normalize) {
    return costMat(x.nrow() % 2, y.nrow()) / fmax(x.nrow(), y.nrow());
  } else {
    return costMat(x.nrow() % 2, y.nrow());
  }
}

//' (Fast) Edit Distance on Real Sequence and Sakoe-Chiba Window
//'
//' Computes the Edit distance on Real Sequence as described by Chen, Özsu
//' and Oria (2005), constraining the possible matches to a maximum index
//' difference of \code{windowSize} as desribed by Sakoe and Chiba (1978).
//' The coding is inspired by the \code{TSdist::EDRDistance()} method, but
//' faster because point-to-point distances computation is integrated into the
//' C++ code.
//'
//' Despite the name, it is not really a distance in the strict sense, as EDR
//' violates the triangular inequality.
//'
//' @section References:
//' Chen, L., Özsu, M. T. & Oria, V. (2005). Robust and fast similarity search
//' for moving object trajectories. In \emph{Proceedings of the 2005 acm sigmod
//' international conference on management of data} (pp. 491–502). ACM.
//'
//' Sakoe, H., & Chiba, S. (1978). Dynamic programming algorithm optimization
//' for spoken word recognition. \emph{IEEE transactions on acoustics, speech,
//' and signal processing, 26}(1), 43-49.
//'
//' @param x 1st numeric vector/time series.
//' @param y 2nd numeric vector/time series.
//' @param epsilon Maximum distance between two time series elements to count a
//' match.
//' @param windowSize The maximum index difference which is considered when
//' matching elements
//' @return The distance as double (not-a-number if matching is not possible
//' as the time series lengths differ more than \code{windowSize}).
//' @family Edit distance functions
//' @export
// [[Rcpp::export]]
double EDRDistSakoeChiba_fast(NumericVector x, NumericVector y, double epsilon, int windowSize) {
  if (fabs(x.size() - y.size()) > windowSize) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  // Create cost matrix; is computed row by row and we only need access to
  // the current plus the previous row, so two rows are sufficient
  // Takes maximum integer - 1, as we add 1 in the minimum calculation below
  // and we don't want an integer overflow
  int maxIntValue = std::numeric_limits<int>::max() - 1;
  IntegerMatrix costMat(2, y.size() + 1);
  // If we are at the right end of the window, value directly above it in
  // imaginary matrix should not be usable -> prevent by filling
  std::fill(costMat.begin(), costMat.end(), maxIntValue);
  // Fill 1st row (only elements within window)
  costMat(0,0) = 0;
  int max1stRowColIndex = fmin(windowSize + 1, y.size());
  for (int j = 1; j <= max1stRowColIndex; j++) {
    costMat(0,j) = j;
  }
  // Find optimal path (considering windowing constraints)
  int minj, maxj, bottom, top;
  for (int i = 1; i <= x.size(); i++) {
    bottom = i % 2; // current row
    top = (i - 1) % 2; // previous row
    minj = fmax(i - windowSize, 1);
    maxj = fmin(i + windowSize, y.size());
    if (minj == 1) { // left window boundary outside matrix
      costMat(bottom, 0) = i;
    } else { // left window boundary inside matrix, we can't use value left of window
      costMat(bottom, minj-1) = maxIntValue;
    }
    for (int j = minj; j <= maxj; j++) {
      costMat(bottom, j) = fmin(fmin(costMat(top, j) + 1, costMat(bottom, j-1) + 1),
              costMat(top, j-1) + (fabs(x[i-1] - y[j-1]) > epsilon));
    }
  }
  return costMat(x.size() % 2, y.size());
}

//' (Fast) Multi-variate Edit Distance on Real Sequence and Sakoe-Chiba Window
//'
//' Computes the Edit distance on Real Sequence as described by Chen, Özsu
//' and Oria (2005), constraining the possible matches to a maximum index
//' difference of \code{windowSize} as desribed by Sakoe and Chiba (1978).
//' The coding is inspired by the \code{TSdist::EDRDistance()} method, but
//' faster because point-to-point distances computation is integrated into the
//' C++ code.
//'
//' Despite the name, it is not really a distance in the strict sense, as EDR
//' violates the triangular inequality.
//'
//' @section References:
//' Chen, L., Özsu, M. T. & Oria, V. (2005). Robust and fast similarity search
//' for moving object trajectories. In \emph{Proceedings of the 2005 acm sigmod
//' international conference on management of data} (pp. 491–502). ACM.
//'
//' Sakoe, H., & Chiba, S. (1978). Dynamic programming algorithm optimization
//' for spoken word recognition. \emph{IEEE transactions on acoustics, speech,
//' and signal processing, 26}(1), 43-49.
//'
//' @param x 1st numeric matrix/multi-variate time series.
//' @param y 2nd numeric matrix/multi-variate time series.
//' @param epsilon Maximum distance between two time series elements to count a
//' match.
//' @param windowSize The maximum index difference which is considered when
//' matching elements
//' @return The distance as double (not-a-number if matching is not possible
//' as the time series lengths differ more than \code{windowSize}).
//' @family Edit distance functions
//' @export
// [[Rcpp::export]]
double EDRDistSakoeChibaMult_fast(NumericMatrix x, NumericMatrix y, double epsilon, int windowSize) {
  if (fabs(x.nrow() - y.nrow()) > windowSize) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  // Create cost matrix; is computed row by row and we only need access to
  // the current plus the previous row, so two rows are sufficient
  // Takes maximum integer - 1, as we add 1 in the minimum calculation below
  // and we don't want an integer overflow
  int maxIntValue = std::numeric_limits<int>::max() - 1;
  IntegerMatrix costMat(2, y.nrow() + 1);
  // If we are at the right end of the window, value directly above it in
  // imaginary matrix should not be usable -> prevent by filling
  std::fill(costMat.begin(), costMat.end(), maxIntValue);
  // Fill 1st row (only elements within window)
  costMat(0,0) = 0;
  int max1stRowColIndex = fmin(windowSize + 1, y.nrow());
  for (int j = 1; j <= max1stRowColIndex; j++) {
    costMat(0,j) = j;
  }
  // Find optimal path (considering windowing constraints)
  int minj, maxj, bottom, top;
  for (int i = 1; i <= x.nrow(); i++) {
    bottom = i % 2; // current row
    top = (i - 1) % 2; // previous row
    minj = fmax(i - windowSize, 1);
    maxj = fmin(i + windowSize, y.nrow());
    if (minj == 1) { // left window boundary outside matrix
      costMat(bottom, 0) = i;
    } else { // left window boundary inside matrix, we can't use value left of window
      costMat(bottom, minj-1) = maxIntValue;
    }
    for (int j = minj; j <= maxj; j++) {
      costMat(bottom, j) = fmin(fmin(costMat(top, j) + 1, costMat(bottom, j-1) + 1),
              costMat(top, j-1) + any(abs(x.row(i-1) - y.row(j-1)) > epsilon).is_true());
    }
  }
  return costMat(x.nrow() % 2, y.nrow());
}

//' (Fast) Edit Distance with Real Penalty
//'
//' Computes the Edit distance with real penalty as described by Chen and Ng
//' (2004). The coding is inspired by the \code{TSdist::ERPDistance()} method,
//' but faster because point-to-point distances computation is integrated into
//' the C++ code.
//'
//' @section References:
//' Chen, L., & Ng, R. (2004, August). On the marriage of lp-norms and edit
//' distance. In \emph{Proceedings of the Thirtieth international conference
//' on Very large data bases-Volume 30} (pp. 792-803). VLDB Endowment.
//'
//' @param x 1st numeric vector/time series.
//' @param y 2nd numeric vector/time series.
//' @param gapValue If an element of one series is not matched to the other
//' series, its distance to the gapValue is computed instead (0 might be a
//' sensible default for standardized series).
//' @return The distance as double.
//' @family Edit distance functions
//' @export
// [[Rcpp::export]]
double ERPDist(NumericVector x, NumericVector y, double gapValue) {
  // Create cost matrix; has one extra column/row (so the index will be one
  // lower when acessing a corresponding element in x and y)
  NumericMatrix costMat(x.size() + 1, y.size() + 1);
  // Fill 1st column and row
  costMat(0,0) = 0.0;
  for (int i = 1; i <= x.size(); i++) {
    costMat(i,0) = costMat(i-1,0) + fabs(x[i-1] - gapValue);
  }
  for (int j = 1; j <= y.size(); j++) {
    costMat(0,j) = costMat(0,j-1) + fabs(y[j-1] - gapValue);
  }
  // Find optimal path
  for (int i = 1; i <= x.size(); i++) {
    for (int j = 1; j <= y.size(); j++) {
      costMat(i,j) = fmin(fmin(costMat(i-1,j) + fabs(x[i-1] - gapValue),
                          costMat(i,j-1) + fabs(y[j-1] - gapValue)),
                          costMat(i-1,j-1) + fabs(x[i-1] - y[j-1]));
    }
  }
  return costMat(x.size(), y.size()); //bottom-right corner (0-indexed)
}

//' (Even Faster) Edit Distance with Real Penalty
//'
//' Faster version of\link{ERPDist} which uses a cyclic access strategy
//' with a smaller cost matrix; inspired by the C implementation of dynamic
//' time warping in \code{dtwclust::dtw_basic()}.
//'
//' @section References:
//' Chen, L., & Ng, R. (2004, August). On the marriage of lp-norms and edit
//' distance. In \emph{Proceedings of the Thirtieth international conference
//' on Very large data bases-Volume 30} (pp. 792-803). VLDB Endowment.
//'
//' @param x 1st numeric vector/time series.
//' @param y 2nd numeric vector/time series.
//' @param gapValue If an element of one series is not matched to the other
//' series, its distance to the gapValue is computed instead (0 might be a
//' sensible default for standardized series).
//' @param normalize Divide by the length of the longer time series (= minimum
//' amount of assignment steps) to account for series of different lengths in
//' your dataset.
//' @return The distance as double.
//' @family Edit distance functions
//' @export
// [[Rcpp::export]]
double ERPDist_fast(NumericVector x, NumericVector y, double gapValue,
                    bool normalize = false) {
  // Create cost matrix; is computed row by row and we only need access to
  // the current plus the previous row, so two rows are sufficient
  NumericMatrix costMat(2, y.size() + 1);
  // Fill 1st row
  costMat(0,0) = 0.0;
  for (int j = 1; j <= y.size(); j++) {
    costMat(0,j) = costMat(0,j-1) + fabs(y[j-1] - gapValue);
  }
  // Find optimal path
  int bottom, top; // we don't shift values when changing rows, just change indices
  for (int i = 1; i <= x.size(); i++) {
    bottom = i % 2; // current row
    top = (i - 1) % 2; // previous row
    costMat(bottom, 0) = costMat(top, 0) + fabs(x[i-1] - gapValue);
    for (int j = 1; j <= y.size(); j++) {
      costMat(bottom, j) = fmin(fmin(costMat(top, j) + fabs(x[i-1] - gapValue),
                                costMat(bottom, j-1) + fabs(y[j-1] - gapValue)),
                                costMat(top, j-1) + fabs(x[i-1] - y[j-1]));
    }
  }
  if (normalize) {
    return costMat(x.size() % 2, y.size()) / fmax(x.size(), y.size());
  } else {
    return costMat(x.size() % 2, y.size());
  }
}

//' (Fast) Multi-variate Edit Distance with Real Penalty
//'
//' Multi-variate version of \code{\link{ERPDist_fast}}. Uses the L1 norm for
//' point-to-point distance computations.
//'
//' @section References:
//' Chen, L., & Ng, R. (2004, August). On the marriage of lp-norms and edit
//' distance. In \emph{Proceedings of the Thirtieth international conference
//' on Very large data bases-Volume 30} (pp. 792-803). VLDB Endowment.
//'
//' @param x 1st numeric matrix/multi-variate time series.
//' @param y 2nd numeric matrix/multi-variate time series.
//' @param gapValue If an element of one series is not matched to the other
//' series, its distance to the gapValue is computed instead (0 might be a
//' sensible default for standardized series).
//' @param normalize Divide by the length of the longer time series (= minimum
//' amount of assignment steps) to account for series of different lengths in
//' your dataset.
//' @return The distance as double.
//' @family Edit distance functions
//' @export
// [[Rcpp::export]]
double ERPDistMult_fast(NumericMatrix x, NumericMatrix y, double gapValue,
                        bool normalize = false) {
  // Create cost matrix; is computed row by row and we only need access to
  // the current plus the previous row, so two rows are sufficient
  NumericMatrix costMat(2, y.nrow() + 1);
  // Fill 1st row
  costMat(0,0) = 0.0;
  for (int j = 1; j <= y.nrow(); j++) {
    costMat(0,j) = costMat(0,j-1) + sum(abs(y.row(j-1) - gapValue));
  }
  // Find optimal path
  int bottom, top; // we don't shift values when changing rows, just change indices
  for (int i = 1; i <= x.nrow(); i++) {
    bottom = i % 2; // current row
    top = (i - 1) % 2; // previous row
    costMat(bottom, 0) = costMat(top, 0) + sum(abs(x.row(i-1) - gapValue));
    for (int j = 1; j <= y.nrow(); j++) {
      costMat(bottom, j) = fmin(fmin(costMat(top, j) + sum(abs(x.row(i-1) - gapValue)),
                                costMat(bottom, j-1) + sum(abs(y.row(j-1) - gapValue))),
                                costMat(top, j-1) + sum(abs(x.row(i-1) - y.row(j-1))));
    }
  }
  if (normalize) {
    return costMat(x.nrow() % 2, y.nrow()) / fmax(x.nrow(), y.nrow());
  } else {
    return costMat(x.nrow() % 2, y.nrow());
  }
}

//' (Fast) Edit Distance with Real Penalty and Sakoe-Chiba Window
//'
//' Computes the Edit distance with real penalty as described by Chen and Ng
//' (2004), constraining the possible matches to a maximum index difference of
//' \code{windowSize} as desribed by Sakoe and Chiba (1978). The coding is
//' inspired by the \code{TSdist::ERPDistance()} method, but faster because
//' point-to-point distances computation is integrated into the C++ code.
//'
//' @section References:
//' Chen, L., & Ng, R. (2004, August). On the marriage of lp-norms and edit
//' distance. In \emph{Proceedings of the Thirtieth international conference
//' on Very large data bases-Volume 30} (pp. 792-803). VLDB Endowment.
//'
//' Sakoe, H., & Chiba, S. (1978). Dynamic programming algorithm optimization
//' for spoken word recognition. \emph{IEEE transactions on acoustics, speech,
//' and signal processing, 26}(1), 43-49.
//'
//' @param x 1st numeric vector/time series.
//' @param y 2nd numeric vector/time series.
//' @param gapValue If an element of one series is not matched to the other
//' series, its distance to the gapValue is computed instead (0 might be a
//' sensible default for standardized series).
//' @param windowSize The maximum index difference which is considered when
//' matching elements.
//' @return The distance as double (not-a-number if matching is not possible
//' as the time series lengths differ more than \code{windowSize}).
//' @family Edit distance functions
//' @export
// [[Rcpp::export]]
double ERPDistSakoeChiba(NumericVector x, NumericVector y, double gapValue, int windowSize) {
  if (fabs(x.size() - y.size()) > windowSize) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  // Create cost matrix, fill with high default values (necessary as elements outside
  // the window should not be considered -> assign high cost); has one extra column/row
  // (so the index will be one lower when acessing a corresponding element in x and y)
  NumericMatrix costMat(x.size() + 1, y.size() + 1);
  std::fill(costMat.begin(), costMat.end(), std::numeric_limits<double>::max());
  // Fill 1st column and row (only elements within window)
  costMat(0,0) = 0.0;
  int max1stColRowIndex = fmin(windowSize + 1, x.size());
  for (int i = 1; i <= max1stColRowIndex; i++) {
    costMat(i,0) = costMat(i-1,0) + fabs(x[i-1] - gapValue);
  }
  int max1stRowColIndex = fmin(windowSize + 1, y.size());
  for (int j = 1; j <= max1stRowColIndex; j++) {
    costMat(0,j) = costMat(0,j-1) + fabs(y[j-1] - gapValue);
  }
  // Find optimal path (considering windowing constraints)
  int minj, maxj;
  for (int i = 1; i <= x.size(); i++) {
    minj = fmax(i - windowSize, 1);
    maxj = fmin(i + windowSize, y.size());
    for (int j = minj; j <= maxj; j++) {
      costMat(i,j) = fmin(fmin(costMat(i-1,j) + fabs(x[i-1] - gapValue),
                          costMat(i,j-1) + fabs(y[j-1] - gapValue)),
                          costMat(i-1,j-1) + fabs(x[i-1] - y[j-1]));
    }
  }
  return costMat(x.size(), y.size()); //bottom-right corner (0-indexed)
}

//' (Even Faster) Edit Distance with Real Penalty and Sakoe-Chiba Window
//'
//' Faster version of \link{ERPDistSakoeChiba} which uses a cyclic access
//' strategy with a smaller cost matrix ; inspired by the C implementation of
//' dynamic time warping in \code{dtwclust::dtw_basic()}.
//'
//' @section References:
//' Chen, L., & Ng, R. (2004, August). On the marriage of lp-norms and edit
//' distance. In \emph{Proceedings of the Thirtieth international conference
//' on Very large data bases-Volume 30} (pp. 792-803). VLDB Endowment.
//'
//' Sakoe, H., & Chiba, S. (1978). Dynamic programming algorithm optimization
//' for spoken word recognition. \emph{IEEE transactions on acoustics, speech,
//' and signal processing, 26}(1), 43-49.
//'
//' @param x 1st numeric vector/time series.
//' @param y 2nd numeric vector/time series.
//' @param gapValue If an element of one series is not matched to the other
//' series, its distance to the gapValue is computed instead (0 might be a
//' sensible default for standardized series).
//' @param windowSize The maximum index difference which is considered when
//' matching elements.
//' @return The distance as double (not-a-number if matching is not possible
//' as the time series lengths differ more than \code{windowSize}).
//' @family Edit distance functions
//' @export
// [[Rcpp::export]]
double ERPDistSakoeChiba_fast(NumericVector x, NumericVector y, double gapValue, int windowSize) {
  if (fabs(x.size() - y.size()) > windowSize) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  // Create cost matrix; is computed row by row and we only need access to
  // the current plus the previous row, so two rows are sufficient
  double maxDoubleValue = std::numeric_limits<double>::max();
  NumericMatrix costMat(2, y.size() + 1);
  // If we are at the right end of the window, value directly above it in
  // imaginary matrix should not be usable -> prevent by filling
  std::fill(costMat.begin(), costMat.end(), maxDoubleValue);
  // Fill 1st row (only elements within window)
  costMat(0,0) = 0.0;
  int max1stRowColIndex = fmin(windowSize + 1, y.size());
  for (int j = 1; j <= max1stRowColIndex; j++) {
    costMat(0,j) = costMat(0,j-1) + fabs(y[j-1] - gapValue);
  }
  // Find optimal path (considering windowing constraints)
  int minj, maxj, bottom, top;
  for (int i = 1; i <= x.size(); i++) {
    bottom = i % 2; // current row
    top = (i - 1) % 2; // previous row
    minj = fmax(i - windowSize, 1);
    maxj = fmin(i + windowSize, y.size());
    if (minj == 1) { // left window boundary outside matrix
      costMat(bottom, 0) = costMat(top, 0) + fabs(x[i-1] - gapValue);
    } else { // left window boundary inside matrix, we can't use value left of window
      costMat(bottom, minj-1) = maxDoubleValue;
    }
    for (int j = minj; j <= maxj; j++) {
      costMat(bottom, j) = fmin(fmin(costMat(top, j) + fabs(x[i-1] - gapValue),
                                costMat(bottom, j-1) + fabs(y[j-1] - gapValue)),
                                costMat(top, j-1) + fabs(x[i-1] - y[j-1]));
    }
  }
  return costMat(x.size() % 2, y.size());
}

//' (Fast) Multi-variate Edit Distance with Real Penalty and Sakoe-Chiba Window
//'
//' Multi-variate version of \link{ERPDistSakoeChiba_fast}. Uses the L1 norm for
//' point-to-point distance computations.
//'
//' @section References:
//' Chen, L., & Ng, R. (2004, August). On the marriage of lp-norms and edit
//' distance. In \emph{Proceedings of the Thirtieth international conference
//' on Very large data bases-Volume 30} (pp. 792-803). VLDB Endowment.
//'
//' Sakoe, H., & Chiba, S. (1978). Dynamic programming algorithm optimization
//' for spoken word recognition. \emph{IEEE transactions on acoustics, speech,
//' and signal processing, 26}(1), 43-49.
//'
//' @param x 1st numeric matrix/multi-variate time series.
//' @param y 2nd numeric vector/multi-variate time series.
//' @param gapValue If an element of one series is not matched to the other
//' series, its distance to the gapValue is computed instead (0 might be a
//' sensible default for standardized series).
//' @param windowSize The maximum index difference which is considered when
//' matching elements.
//' @return The distance as double (not-a-number if matching is not possible
//' as the time series lengths differ more than \code{windowSize}).
//' @family Edit distance functions
//' @export
// [[Rcpp::export]]
double ERPDistSakoeChibaMult_fast(NumericMatrix x, NumericMatrix y, double gapValue, int windowSize) {
  if (fabs(x.nrow() - y.nrow()) > windowSize) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  // Create cost matrix; is computed row by row and we only need access to
  // the current plus the previous row, so two rows are sufficient
  double maxDoubleValue = std::numeric_limits<double>::max();
  NumericMatrix costMat(2, y.nrow() + 1);
  // If we are at the right end of the window, value directly above it in
  // imaginary matrix should not be usable -> prevent by filling
  std::fill(costMat.begin(), costMat.end(), maxDoubleValue);
  // Fill 1st row (only elements within window)
  costMat(0,0) = 0.0;
  int max1stRowColIndex = fmin(windowSize + 1, y.nrow());
  for (int j = 1; j <= max1stRowColIndex; j++) {
    costMat(0,j) = costMat(0,j-1) + sum(abs(y.row(j-1) - gapValue));
  }
  // Find optimal path (considering windowing constraints)
  int minj, maxj, bottom, top;
  for (int i = 1; i <= x.nrow(); i++) {
    bottom = i % 2; // current row
    top = (i - 1) % 2; // previous row
    minj = fmax(i - windowSize, 1);
    maxj = fmin(i + windowSize, y.nrow());
    if (minj == 1) { // left window boundary outside matrix
      costMat(bottom, 0) = costMat(top, 0) + sum(abs(x.row(i-1) - gapValue));
    } else { // left window boundary inside matrix, we can't use value left of window
      costMat(bottom, minj-1) = maxDoubleValue;
    }
    for (int j = minj; j <= maxj; j++) {
      costMat(bottom, j) = fmin(fmin(costMat(top, j) + sum(abs(x.row(i-1) - gapValue)),
                                costMat(bottom, j-1) + sum(abs(y.row(j-1) - gapValue))),
                                costMat(top, j-1) + sum(abs(x.row(i-1) - y.row(j-1))));
    }
  }
  return costMat(x.nrow() % 2, y.nrow());
}

//' (Fast) Dynamic Time Warping Dissimilarity
//'
//' Fast version of univariate dynamic time warping (unconstrained, symmetric1
//' step pattern) which uses a cyclic access strategy with a smaller cost matrix;
//' inspired by the C implementation of dynamic time warping in
//' \code{dtwclust::dtw_basic()}, but cuts even more overhead.
//'
//' Be aware that it is not really a distance in the strict sense, as DTW
//' violates the triangle inequality.
//'
//' @section References:
//' Berndt, D. J. & Clifford, J. (1994). Using dynamic time warping to find
//' patterns in time series. In \emph{Proceedings of the 3rd international
//' conference on knowledge discovery and data mining} (pp. 359–370). AAAI Press.
//'
//' @param x 1st numeric vector/time series.
//' @param y 2nd numeric vector/time series.
//' @param cid Should the distance be made "complexity-invariant"
//' (\code{\link{l2CompCorFactor_fast}})?
//' @param cortK Should the temporal behavior (correlation) of the time series'
//' diff vectors be considered (\code{\link{cortFactor_fast}})? A factor smaller
//' than 0 means no, higher factors will be used as parameter \code{k} in the
//' temporal correlation scaling function.
//' @param normalize Divide by the length of the longer time series (= minimum
//' amount of assignment steps) to account for series of different lengths in
//' your dataset.
//' @return The distance as double.
//' @family DTW functions
//' @export
// [[Rcpp::export]]
double DTWDist_fast(NumericVector x, NumericVector y, bool cid = false,
                    double cortK = -1, bool normalize = false) {
  // Create cost matrix; is computed row by row and we only need access to
  // the current plus the previous row, so two rows are sufficient
  NumericMatrix costMat(2, y.size() + 1);
  // Fill 1st row (also second, but this is overwritten anyway)
  double maxDoubleValue = std::numeric_limits<double>::max();
  std::fill(costMat.begin(), costMat.end(), maxDoubleValue);
  costMat(0,0) = 0.0;
  // Find optimal path
  int bottom, top; // we don't shift values when changing rows, just change indices
  for (int i = 1; i <= x.size(); i++) {
    bottom = i % 2; // current row
    top = (i - 1) % 2; // previous row
    if (i == 2) { //only for 1st iteration interesting as default
      costMat(0,0) = maxDoubleValue;
    }
    for (int j = 1; j <= y.size(); j++) {
      costMat(bottom, j) = fabs(x[i-1] - y[j-1]) +
        fmin(fmin(costMat(top, j), costMat(bottom, j-1)), costMat(top, j-1));
    }
  }
  double result = costMat(x.size() % 2, y.size());
  if (normalize) {
    result /= fmax(x.size(), y.size());
  }
  if (cid) {
    result *= l2CompCorFactor_fast(x,y);
  }
  if (cortK >= 0) {
    result *= cortFactor_fast(x,y);
  }
  return result;
}

//' (Fast) Multi-variate Dynamic Time Warping Dissimilarity
//'
//' Fast version of multi-variate dynamic time warping (unconstrained, symmetric1
//' step pattern, L2 distance for point-to-point comparisons) which uses a cyclic
//' access strategy with a smaller cost matrix; inspired by the C implementation of
//' dynamic time warping in \code{dtwclust::dtw_basic()}, but cuts even more overhead.
//'
//' Be aware that it is not really a distance in the strict sense, as DTW
//' violates the triangle inequality.
//'
//' @section References:
//'
//' Berndt, D. J. & Clifford, J. (1994). Using dynamic time warping to find
//' patterns in time series. In \emph{Proceedings of the 3rd international
//' conference on knowledge discovery and data mining} (pp. 359–370). AAAI Press.
//'
//' @param x 1st numeric matrix/multi-variate time series.
//' @param y 2nd numeric matrix/multi-variate time series.
//' @param cid Should the distance be made "complexity-invariant"
//' (\code{\link{l2CompCorFactorMult_fast}})?
//' @param cortK Should the temporal behavior (correlation) of the time series'
//' diff vectors be considered (\code{\link{cortFactorMult_fast}})? A factor
//' smaller than 0 means no, higher factors will be used as parameter \code{k}
//' in the temporal correlation scaling function.
//' @param normalize Divide by the length of the longer time series (= minimum
//' amount of assignment steps) to account for series of different lengths in
//' your dataset.
//' @return The distance as double.
//' @family DTW functions
//' @export
// [[Rcpp::export]]
double DTWDistMult_fast(NumericMatrix x, NumericMatrix y, bool cid = false,
                        double cortK = -1, bool normalize = false) {
  // Create cost matrix; is computed row by row and we only need access to
  // the current plus the previous row, so two rows are sufficient
  NumericMatrix costMat(2, y.nrow() + 1);
  // Fill 1st row (also second, but this is overwritten anyway)
  double maxDoubleValue = std::numeric_limits<double>::max();
  std::fill(costMat.begin(), costMat.end(), maxDoubleValue);
  costMat(0,0) = 0.0;
  // Find optimal path
  int bottom, top; // we don't shift values when changing rows, just change indices
  for (int i = 1; i <= x.nrow(); i++) {
    bottom = i % 2; // current row
    top = (i - 1) % 2; // previous row
    if (i == 2) { //only for 1st iteration interesting as default
      costMat(0,0) = maxDoubleValue;
    }
    for (int j = 1; j <= y.nrow(); j++) {
      costMat(bottom, j) = sqrt(sum(pow(x.row(i-1) - y.row(j-1), 2))) +
        fmin(fmin(costMat(top, j), costMat(bottom, j-1)), costMat(top, j-1));
    }
  }
  double result = costMat(x.nrow() % 2, y.nrow());
  if (normalize) {
    result /= fmax(x.nrow(), y.nrow());
  }
  if (cid) {
    result *= l2CompCorFactorMult_fast(x,y);
  }
  if (cortK >= 0) {
    result *= cortFactorMult_fast(x,y);
  }
  return result;
}

//' (Fast) Dynamic Time Warping Dissimilarity with a Sakoe-Chiba Window
//'
//' Fast version of univariate dynamic time warping (Sakoe-Chiba window as
//' constraint, symmetric1 step pattern) which uses a cyclic access strategy
//' with a smaller cost matrix; inspired by the C implementation of dynamic time
//' warping in \code{dtwclust::dtw_basic()}, but cuts even more overhead.
//'
//' Be aware that it is not really a distance in the strict sense, as DTW
//' violates the triangle inequality.
//'
//' @section References:
//' Berndt, D. J. & Clifford, J. (1994). Using dynamic time warping to find
//' patterns in time series. In \emph{Proceedings of the 3rd international
//' conference on knowledge discovery and data mining} (pp. 359–370). AAAI Press.
//'
//' Sakoe, H., & Chiba, S. (1978). Dynamic programming algorithm optimization
//' for spoken word recognition. \emph{IEEE transactions on acoustics, speech,
//' and signal processing, 26}(1), 43-49.
//'
//' @param x 1st numeric vector/time series.
//' @param y 2nd numeric vector/time series.
//' @param windowSize The maximum index difference which is considered when
//' matching elements. If greater/equal one, interpreted as absolute value.
//' If smaller then one, interpreted as fraction of the length of the longer
//' time series.
//' @return The distance as double (not-a-number if matching is not possible
//' as the time series lengths differ more than \code{windowSize}).
//' @family DTW functions
//' @export
// [[Rcpp::export]]
double DTWDistSakoeChiba_fast(NumericVector x, NumericVector y, double windowSize) {
  int intWindowSize;
  if (windowSize >= 1) {
    intWindowSize = std::floor(windowSize);
  } else {
    intWindowSize = std::floor(std::max(x.size(), y.size()) * windowSize);
  }
  if (fabs(x.size() - y.size()) > intWindowSize) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  // Create cost matrix; is computed row by row and we only need access to
  // the current plus the previous row, so two rows are sufficient
  NumericMatrix costMat(2, y.size() + 1);
  // Fill 1st row (also second, but this is overwritten anyway)
  double maxDoubleValue = std::numeric_limits<double>::max();
  std::fill(costMat.begin(), costMat.end(), std::numeric_limits<double>::max());
  costMat(0,0) = 0.0;
  // Find optimal path (considering windowing constraints)
  int minj, maxj, bottom, top;
  for (int i = 1; i <= x.size(); i++) {
    bottom = i % 2; // current row
    top = (i - 1) % 2; // previous row
    minj = fmax(i - intWindowSize, 1);
    maxj = fmin(i + intWindowSize, y.size());
    if (i == 2) { //only for 1st iteration interesting as default
      costMat(0,0) = maxDoubleValue;
    }
    if (minj > 1) { // left window boundary inside matrix, we can't use value left of window
      costMat(bottom, minj-1) = maxDoubleValue;
    }
    for (int j = minj; j <= maxj; j++) {
      costMat(bottom, j) = fabs(x[i-1] - y[j-1]) +
        fmin(fmin(costMat(top, j), costMat(bottom, j-1)), costMat(top, j-1));
    }
  }
  return costMat(x.size() % 2, y.size());
}

//' (Fast) Multi-variate Dynamic Time Warping Dissimilarity with a Sakoe-Chiba Window
//'
//' Fast version of multivariate dynamic time warping (Sakoe-Chiba window as
//' constraint, symmetric1 step pattern, L2 distance for point-to-point comparisons)
//' which uses a cyclic access strategy with a smaller cost matrix; inspired by
//' the C implementation of dynamic time warping in \code{dtwclust::dtw_basic()},
//' but cuts even more overhead.
//'
//' Be aware that it is not really a distance in the strict sense, as DTW
//' violates the triangle inequality.
//'
//' @section References:
//' Berndt, D. J. & Clifford, J. (1994). Using dynamic time warping to find
//' patterns in time series. In \emph{Proceedings of the 3rd international
//' conference on knowledge discovery and data mining} (pp. 359–370). AAAI Press.
//'
//' Sakoe, H., & Chiba, S. (1978). Dynamic programming algorithm optimization
//' for spoken word recognition. \emph{IEEE transactions on acoustics, speech,
//' and signal processing, 26}(1), 43-49.
//'
//' @param x 1st numeric matrix/multi-variate time series.
//' @param y 2nd numeric matrix/multi-variate time series.
//' @param windowSize The maximum index difference which is considered when
//' matching elements. If greater/equal one, interpreted as absolute value.
//' If smaller then one, interpreted as fraction of the length of the longer
//' time series.
//' @return The distance as double (not-a-number if matching is not possible
//' as the time series lengths differ more than \code{windowSize}).
//' @family DTW functions
//' @export
// [[Rcpp::export]]
double DTWDistSakoeChibaMult_fast(NumericMatrix x, NumericMatrix y, double windowSize) {
  int intWindowSize;
  if (windowSize >= 1) {
    intWindowSize = std::floor(windowSize);
  } else {
    intWindowSize = std::floor(std::max(x.size(), y.size()) * windowSize);
  }
  if (fabs(x.nrow() - y.nrow()) > intWindowSize) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  // Create cost matrix; is computed row by row and we only need access to
  // the current plus the previous row, so two rows are sufficient
  NumericMatrix costMat(2, y.nrow() + 1);
  // Fill 1st row (also second, but this is overwritten anyway)
  double maxDoubleValue = std::numeric_limits<double>::max();
  std::fill(costMat.begin(), costMat.end(), std::numeric_limits<double>::max());
  costMat(0,0) = 0.0;
  // Find optimal path (considering windowing constraints)
  int minj, maxj, bottom, top;
  for (int i = 1; i <= x.nrow(); i++) {
    bottom = i % 2; // current row
    top = (i - 1) % 2; // previous row
    minj = fmax(i - intWindowSize, 1);
    maxj = fmin(i + intWindowSize, y.nrow());
    if (i == 2) { //only for 1st iteration interesting as default
      costMat(0,0) = maxDoubleValue;
    }
    if (minj > 1) { // left window boundary inside matrix, we can't use value left of window
      costMat(bottom, minj-1) = maxDoubleValue;
    }
    for (int j = minj; j <= maxj; j++) {
      costMat(bottom, j) = sqrt(sum(pow(x.row(i-1) - y.row(j-1), 2))) +
        fmin(fmin(costMat(top, j), costMat(bottom, j-1)), costMat(top, j-1));
    }
  }
  return costMat(x.nrow() % 2, y.nrow());
}

//' (Fast) L1 Distance
//'
//' Computes the standard Manhattan distance with a fast C++ implementation.
//'
//' @param x 1st numeric vector/time series.
//' @param y 2nd numeric vector/time series.
//' @return The distance as double.
//' @family L_p distances
//' @export
// [[Rcpp::export]]
double l1Dist_fast(NumericVector x, NumericVector y) {
  double result = 0;
  int n = x.size();
  for (int i = 0; i < n; i++) {
    result += fabs(x[i]-y[i]);
  }
  return result;
}

//' (Fast) Multi-variate L1 Distance
//'
//' Computes the Manhattan distance between multi-variate time series (according
//' to Kotsakos, Trajcevski, Gunopulos and Aggarwal (2014)) with a fast C++
//' implementation.
//'
//' @section References:
//'
//' Kotsakos, D., Trajcevski, G., Gunopulos, D. & Aggarwal, C. C. (2014). Time-series
//' data clustering. In C. C. Aggarwal & C. K. Reddy (Eds.), \emph{Data clustering :
//' Algorithms and applications} (pp. 357–380). Chapman & Hall/CRC data mining and
//' knowledge discovery series. Boca Raton: CRC Press.
//'
//' @param x 1st numeric matrix/multi-variate time series.
//' @param y 2nd numeric matrix/multi-variate time series.
//' @return The distance as double.
//' @family L_p distances
//' @export
// [[Rcpp::export]]
double l1DistMult_fast(NumericMatrix x, NumericMatrix y) {
  return sum(abs(x - y));
}

//' (Fast) L2 Distance
//'
//' Computes the standard Euclidean distance with a fast C++ implementation.
//'
//' @param x 1st numeric vector/time series.
//' @param y 2nd numeric vector/time series.
//' @param cid Should the distance be made complexity invariant
//' (\code{\link{l2CompCorFactor_fast}})?
//' @param cortK Should the temporal behavior (correlation) of the time series'
//' diff vectors be considered (\code{\link{cortFactor_fast}})? A factor smaller
//' than 0 means no, higher factors will be used as parameter \code{k} in the
//' temporal correlation scaling function.
//' @return The distance as double.
//' @family L_p distances
//' @export
// [[Rcpp::export]]
double l2Dist_fast(NumericVector x, NumericVector y, bool cid = false, double cortK = -1) {
  double result = 0, absDiff;
  int n = x.size();
  for (int i = 0; i < n; i++) {
    absDiff = fabs(x[i]-y[i]);
    result += absDiff*absDiff;
  }
  result = sqrt(result);
  if (cid) {
    return result *= l2CompCorFactor_fast(x,y);
  }
  if (cortK >= 0) {
    result *= cortFactor_fast(x,y,cortK);
  }
  return result;
}

//' (Fast) Multi-variate L2 Distance
//'
//' Computes the standard Euclidean distance between multi-variate time series
//' (according to Kotsakos, Trajcevski, Gunopulos and Aggarwal (2014)) with a
//' fast C++ implementation.
//'
//' @section References:
//'
//' Kotsakos, D., Trajcevski, G., Gunopulos, D. & Aggarwal, C. C. (2014). Time-series
//' data clustering. In C. C. Aggarwal & C. K. Reddy (Eds.), \emph{Data clustering :
//' Algorithms and applications} (pp. 357–380). Chapman & Hall/CRC data mining and
//' knowledge discovery series. Boca Raton: CRC Press.
//'
//' @param x 1st numeric matrix/multi-variate time series.
//' @param y 2nd numeric matrix/multi-variate time series.
//' @param cid Should the distance be made complexity invariant
//' (\code{\link{l2CompCorFactorMult_fast}})?
//' @param cortK Should the temporal behavior (correlation) of the time series'
//' diff vectors be considered (\code{\link{cortFactorMult_fast}})? A factor smaller
//' than 0 means no, higher factors will be used as parameter \code{k} in the
//' temporal correlation scaling function.
//' @return The distance as double.
//' @family L_p distances
//' @export
// [[Rcpp::export]]
double l2DistMult_fast(NumericMatrix x, NumericMatrix y, bool cid = false, double cortK = -1) {
  double result = sqrt(sum((x-y) * (x-y)));
  if (cid) {
    return result *= l2CompCorFactorMult_fast(x,y);
  }
  if (cortK >= 0) {
    result *= cortFactorMult_fast(x,y,cortK);
  }
  return result;
}

//' (Fast) Chebyshev Distance
//'
//' Computes the standard Chebyshev distance (maximum metric) with a fast C++
//' implementation.
//'
//' @param x 1st numeric vector/time series.
//' @param y 2nd numeric vector/time series.
//' @return The distance as double.
//' @family L_p distances
//' @export
// [[Rcpp::export]]
double lmaxDist_fast(NumericVector x, NumericVector y) {
  double result = 0;
  int n = x.size();
  for (int i = 0; i < n; i++) {
    result = fmax(result, fabs(x[i]-y[i]));
  }
  return result;
}

//' (Fast) Multi-variate Chebyshev Distance
//'
//' Computes the standard Chebyshev distance (maximum metric) between
//' multi-variate time series (according toKotsakos, Trajcevski, Gunopulos and
//' Aggarwal (2014)) with a fast C++ implementation.
//'
//' @section References:
//'
//' Kotsakos, D., Trajcevski, G., Gunopulos, D. & Aggarwal, C. C. (2014). Time-series
//' data clustering. In C. C. Aggarwal & C. K. Reddy (Eds.), \emph{Data clustering :
//' Algorithms and applications} (pp. 357–380). Chapman & Hall/CRC data mining and
//' knowledge discovery series. Boca Raton: CRC Press.
//'
//' @param x 1st numeric matrix/multi-variate time series.
//' @param y 2nd numeric matrix/multi-variate time series.
//' @return The distance as double.
//' @family L_p distances
//' @export
// [[Rcpp::export]]
double lmaxDistMult_fast(NumericMatrix x, NumericMatrix y) {
  return max(abs(x - y));
}
