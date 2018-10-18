#include <Rcpp.h>
using namespace Rcpp;

//' Entropy for Clusterings
//'
//' Calculates the Shannon entropy for a cluster assignment vector. A value of
//' 0 means that all elements are in one cluster, higher values indicate a more
//' even distribution of objects in the clusters. The value can be normalized
//' to [0,1] such that 1 means the same number of objects in each cluster.
//'
//' @param assignments Integer vector of cluster assignments containing only values
//' from 1 to k with k = number of clusters (code depends on this!).
//' @param normalize Should the entropy be normalized to [0,1]?
//' @return The entropy as double.
//' @export
// [[Rcpp::export]]
double clusterEntropy_fast(IntegerVector assignments, bool normalize = false) {
  NumericVector table = NumericVector(unique(assignments).size());
  for (int i = 0; i < assignments.size(); i++) {
    table[assignments[i] - 1]++;
  }
  table = table / assignments.size();
  double sum = 0;
  for (int i = 0; i < table.size(); i++) {
    if (table[i] != 0) {
      sum -= table[i] * log2(table[i]);
    }
  }
  if (normalize) {
    if (table.size() == 0 || table.size() == 1) {
      return 0;
    } else {
      return sum / log2(table.size());
    }
  } else {
    return sum;
  }
}

//' Generalized Davies-Bouldin Index
//'
//' Calculates a generalized version of the Davies-Bouldin Index for internal
//' cluster validation. The index is expressed by a ratio of cluster compactness
//' and cluster separation, summed and averaged over all clusters. This generalized
//' method does not define the concrete way to compute distances within clusters
//' and between clusters, but simply takes these distances as input to compute the
//' Davies-Bouldin Index. Lower values indicate better clustering quality.
//'
//' @section References:
//'
//' Davies, D. L. & Bouldin, D. W. (1979). A cluster separation measure. \emph{IEEE
//' transactions on pattern analysis and machine intelligence, 1}(2), 224–227.
//'
//' @param interClusterDistances A matrix representing the distances between
//' clusters (separation), with the number of rows/columns equal to the number
//' of clusters.
//' @param intraClusterDistances A vector representing the distances within a
//' cluster (compactness), with its length equal to the number of clusters.
//' @return The Generalized Davies-Bouldin Index. Could be Inf or NaN if there
//' are clusters with a distance of zero between them (which is a bad clustering
//' result).
//' @family Internal Cluster Validity Indices
//' @export
// [[Rcpp::export]]
double generalizedDB_fast(NumericMatrix interClusterDistances, NumericVector intraClusterDistances) {
  double result = 0, maxForCluster;
  for (int i = 0; i < interClusterDistances.nrow(); i++) {
    maxForCluster = 0;
    for (int j = 0; j < interClusterDistances.ncol(); j++) {
      if (i != j) {
        if (interClusterDistances(i,j) != 0) {
          // for each cluster: take "worst" (aka max, DB is to be minimized) intra/inter ratio
          maxForCluster = fmax(maxForCluster, (intraClusterDistances[i] + intraClusterDistances[j]) /
            interClusterDistances(i,j));
        } else if (intraClusterDistances[i] == 0 && intraClusterDistances[j] == 0) {
          return std::numeric_limits<double>::quiet_NaN();
        } else {
          return std::numeric_limits<double>::infinity();
        }
      }
    }
    result += maxForCluster;
  }
  return result / intraClusterDistances.size();
}

//' Inverted Generalized Davies-Bouldin Index
//'
//' Calculates a generalized version of the Davies-Bouldin Index, similar to
//' \code{\link{generalizedDB_fast}}. The only difference is that the separation
//' measure is divided by the compactness measure (inverted compared to original
//' index), so high values are desirable.
//'
//' @section References:
//'
//' Davies, D. L. & Bouldin, D. W. (1979). A cluster separation measure. \emph{IEEE
//' transactions on pattern analysis and machine intelligence, 1}(2), 224–227.
//'
//' @param interClusterDistances A matrix representing the distances between
//' clusters (separation), with the number of rows/columns equal to the number
//' of clusters.
//' @param intraClusterDistances A vector representing the distances within a
//' cluster (compactness), with its length equal to the number of clusters.
//' @return The Inverted Generalized Davies-Bouldin Index. Could be Inf or NaN
//' if there are clusters with a distance of zero between them (which is a bad
//' clustering result) or infinity between them (which might indicate an error
//' in the computation of your dissimilarity matrix).
//' @family Internal Cluster Validity Indices
//' @export
// [[Rcpp::export]]
double iGeneralizedDB_fast(NumericMatrix interClusterDistances, NumericVector intraClusterDistances) {
  double result = 0, minForCluster;
  for (int i = 0; i < interClusterDistances.nrow(); i++) {
    minForCluster = std::numeric_limits<double>::infinity();
    for (int j = 0; j < interClusterDistances.ncol(); j++) {
      if (i != j) {
        if (intraClusterDistances[i] + intraClusterDistances[j] != 0) {
          // for each cluster: take "worst" (aka min, inv DB index is to be maximized) inter/intra ratio
          minForCluster = fmin(minForCluster, interClusterDistances(i,j) /
            (intraClusterDistances[i] + intraClusterDistances[j]));
        } else if (interClusterDistances(i,j) == 0) {
          return std::numeric_limits<double>::quiet_NaN();
        } else {
          return std::numeric_limits<double>::infinity();
        }
      }
    }
    if (!std::isfinite(minForCluster)) {
      return std::numeric_limits<double>::infinity();
    } else {
      result += minForCluster;
    }
  }
  return result / intraClusterDistances.size();
}

//' Generalized Dunn Index
//'
//' Calculates a generalized version of the Dunn Index, allowing an arbitrary
//' measure of cluster separation (which goes to the numerator and is minimized)
//' and an arbitrary measure of cluster compactness (which goes to the denominator
//' and is maximized over all clusters). Dunn used the highly outlier-prone
//' single linkage measure for separation and complete linkage for compactness.
//' Higher values indicate better clustering quality.
//'
//' @section References:
//'
//' Bezdek, J. C. & Pal, N. R. (1998). Some new indexes of cluster validity.
//' \emph{IEEE Transactions on Systems, Man, and Cybernetics, Part B (Cybernetics),
//' 28}(3), 301–315.
//'
//' Dunn, J. C. (1973). A fuzzy relative of the isodata process and its use in
//' detecting compact well-separated clusters. \emph{Journal of Cybernetics, 3}(3),
//' 32–57.
//'
//' @param interClusterDistances A symmetric matrix representing the distances
//' between clusters (separation), with the number of rows/columns equal to the
//' number of clusters.
//' @param intraClusterDistances A vector representing the distances within a
//' cluster (compactness), with its length equal to the number of clusters.
//' @return The Generalized Dunn Index. Could be Inf or NaN if there are clusters
//' with a distance of zero between them (which is a bad clustering result).
//' @family Internal Cluster Validity Indices
//' @export
// [[Rcpp::export]]
double generalizedDunn_fast(NumericMatrix interClusterDistances, NumericVector intraClusterDistances) {
  double minInterClusterDistance = std::numeric_limits<double>::max(),
    maxIntraClusterDistance = max(intraClusterDistances);
  for (int i = 0; i < interClusterDistances.nrow(); i++) {
    for (int j = i + 1; j < interClusterDistances.ncol(); j++) {
      minInterClusterDistance = fmin(minInterClusterDistance, interClusterDistances(i,j));
    }
  }
  if (maxIntraClusterDistance == 0) {
    if (minInterClusterDistance == 0) {
      return std::numeric_limits<double>::quiet_NaN();
    } else {
      return std::numeric_limits<double>::infinity();
    }
  } else {
    return minInterClusterDistance / maxIntraClusterDistance;
  }
}

//' Conditional Entropy to Compare Clusterings
//'
//' Calculates the conditional Shannon entropy to compare two cluster assignment
//' vectors (external cluster validation). It is a value greater or equal 0,
//' lower values indicating more similarity (purer clusters).  Optionally, the
//' index can be normalized to [0,1] and we take 1 - normalized entropy to get
//' a uniformity measure where high values are good.
//'
//' Be aware that this measure is asymmetric (classes are conditioned on/ analyzed
//' in) clusters and can still be high if the classes of the ground truth are
//' split up into multiple (but pure) clusters. Wu, Xiong and Chen (2009) propose
//' to use the symmetric variation of information (\code{\link{VI_fast}}) instead,
//' which is also based on entropy.
//'
//' We use the base 2 logarithm for calculating entropy.
//'
//' @section References:
//'
//' Wu, J., Xiong, H. & Chen, J. (2009). Adapting the right measures for k-means
//' clustering. In \emph{Proceedings of the 15th acm sigkdd international conference
//' on knowledge discovery and data mining} (pp. 877-886). ACM.
//'
//' @param assignments Integer vector of cluster assignments containing only values
//' from 1 to k with k = number of clusters (code depends on this!).
//' @param groundTruth Integer vector of true class labels containing only values
//' from 1 to k' with k' = number of classes.
//' @param normalizeAndInvert Should the entropy be normalized to [0,1] and inverted
//' such that high values indicate similar clusterings?
//' @return The conditional entropy as double in [0, k'] (without normalization) or
//' a uniformity measure in [0,1] (with normalization).
//' @family External Cluster Validity Indices
//' @export
// [[Rcpp::export]]
double conditionalEntropy_fast(IntegerVector assignments, IntegerVector groundTruth,
               bool normalizeAndInvert = false) {
  NumericVector assignmentsTable = NumericVector(unique(assignments).size());
  NumericMatrix sharedTable = NumericMatrix(assignmentsTable.size(), unique(groundTruth).size());
  if (sharedTable.ncol() == 1) { // only one class, clusters always pure
    if (normalizeAndInvert) {
      return 1;
    } else {
      return 0;
    }
  }
  for (int i = 0; i < assignments.size(); i++) {
    assignmentsTable[assignments[i] - 1]++;
    sharedTable(assignments[i] - 1, groundTruth[i] - 1)++;
  }
  assignmentsTable = assignmentsTable / assignments.size();
  sharedTable = sharedTable / assignments.size();
  double entropy = 0, clusterEntropy = 0;
  for (int i = 0; i < sharedTable.nrow(); i++) {
    clusterEntropy = 0;
    for (int j = 0; j < sharedTable.ncol(); j++) {
      if (sharedTable(i,j) != 0) {
        clusterEntropy += sharedTable(i,j) / assignmentsTable[i] *
          log2(sharedTable(i,j) / assignmentsTable[i]);
      }
    }
    entropy -= assignmentsTable[i] * clusterEntropy;
  }
  if (normalizeAndInvert) {
    return 1 - entropy / log2(sharedTable.ncol());
  } else {
    return entropy;
  }
}

//' Statistics for External CVIs based on Pairwise Comparison
//'
//' Calculates the summary statistics [m, m1, m2, M] which can be used to compute
//' multiple normalized external CVIs based on the formulas of Wu, Xiong and Chen (2009).
//'
//' @section References:
//'
//' Wu, J., Xiong, H. & Chen, J. (2009). Adapting the right measures for k-means
//' clustering. In \emph{Proceedings of the 15th acm sigkdd international conference
//' on knowledge discovery and data mining} (pp. 877-886). ACM.
//'
//' @param assignments1 Integer vector of cluster assignments containing only values
//' from 1 to k with k = number of clusters (code depends on this!).
//' @param assignments2 Integer vector of cluster assignments containing only values
//' from 1 to k with k = number of clusters.
//' @return Vector with the four components called \code{m, m1, m2, M} by Wu,
//' Xiong and Chen (2009).
//' @family External Cluster Validity Indices
//' @export
// [[Rcpp::export]]
NumericVector pairCVIParameters_fast(IntegerVector assignments1, IntegerVector assignments2) {
  IntegerVector table1 = IntegerVector(unique(assignments1).size());
  IntegerVector table2 = IntegerVector(unique(assignments2).size());
  IntegerMatrix table12 = IntegerMatrix(table1.size(), table2.size());
  for (int i = 0; i < assignments1.size(); i++) {
    table1[assignments1[i] - 1]++;
    table2[assignments2[i] - 1]++;
    table12(assignments1[i] - 1, assignments2[i] - 1)++;
  }
  NumericVector result = NumericVector(4); // [m, m1, m2, M]
  result[3] = 0.5 * assignments1.size() * (assignments1.size() - 1);
  for (int i = 0; i < table1.size(); i++) {
    result[1] += table1[i] * (table1[i] - 1);
  }
  result[1] /= 2;
  for (int j = 0; j < table2.size(); j++) {
    result[2] += table2[j] * (table2[j] - 1);
  }
  result[2] /= 2;
  for (int i = 0; i < table12.nrow(); i++) {
    for (int j = 0; j < table12.ncol(); j++) {
      result[0] += table12(i,j) * (table12(i,j) - 1);
    }
  }
  result[0] /= 2;
  return result;
}

//' Rand Index
//'
//' Calculates the Rand Index of Rand (1971) to compare two cluster assignment
//' vectors (external cluster validation). It is a value in (0,1], higher values
//' indicating more similarity. The index can be corrected for similarity by
//' chance as proposed by Hubert and Arabie (1985), then also possibly yielding
//' negative results (if the similarity is worse than random assignment) while
//' the maximum is still 1 (and values are usually positive).
//'
//' @section References:
//'
//' Hubert, L. & Arabie, P. (1985). Comparing partitions. \emph{Journal of
//' classification,2}(1), 193-218.
//'
//' Rand, W. M. (1971). Objective criteria for the evaluation of clustering methods.
//' \emph{Journal of the American Statistical association, 66}(336), 846-850.
//'
//' Wu, J., Xiong, H. & Chen, J. (2009). Adapting the right measures for k-means
//' clustering. In \emph{Proceedings of the 15th acm sigkdd international conference
//' on knowledge discovery and data mining} (pp. 877-886). ACM.
//'
//' @param pairCVIParams Output of \code{link{pairCVIParameters}} which has to be
//' called with the two cluster assignment vector to be compared.
//' @param normalize Should the Rand Index be corrected for chance? (Adjusted Rand
//' Index as proposed by Hubert and Arabie (1985))
//' @return The (Adjusted) Rand Index as double (at most 1 for identical clusterings,
//' normal Rand Index greater than zero, adjusted one can also be negative).
//' @family External Cluster Validity Indices
//' @export
// [[Rcpp::export]]
double randIndex_fast(NumericVector pairCVIParams, bool normalize = false) {
  if (normalize) {
    return (pairCVIParams[0] - pairCVIParams[1]*pairCVIParams[2]/pairCVIParams[3]) /
      (0.5*pairCVIParams[1] + 0.5*pairCVIParams[2] - pairCVIParams[1]*pairCVIParams[2]/pairCVIParams[3]);
  } else {
    return (pairCVIParams[3] - pairCVIParams[1] - pairCVIParams[2] + 2*pairCVIParams[0]) /
      pairCVIParams[3];
  }
}

//' Fowlkes-Mallows Index
//'
//' Calculates the index of Fowlkes and Mallows (1983) to compare two cluster
//' assignment vectors (external cluster validation). It is a value in [0,1],
//' the geometric mean of precision and recall, higher values indicating more
//' similarity. The index can be corrected for similarity by chance as proposed
//' by Wu, Xiong and Chen (2009), then also possibly yielding negative results
//' (if the similarity is worse than random assignment) while the maximum is still
//' 1 (and values are usually positive).
//'
//' @section References:
//'
//' Fowlkes, E. B. & Mallows, C. L. (1983). A method for comparing two hierarchical
//' clusterings. \emph{Journal of the American statistical association, 78}(383),
//' 553-569.
//'
//' Wu, J., Xiong, H. & Chen, J. (2009). Adapting the right measures for k-means
//' clustering. In \emph{Proceedings of the 15th acm sigkdd international conference
//' on knowledge discovery and data mining} (pp. 877-886). ACM.
//'
//' @param pairCVIParams Output of \code{link{pairCVIParameters}} which has to be
//' called with the two cluster assignment vector to be compared.
//' @param normalize Should the Fowlkes-Mallows Index be corrected for chance?
//' @return The Fowlkes-Mallows Index as double (at most 1 for identical clusterings,
//' normal Fowlkes-Mallows greater than zero, normalized one can also be negative).
//' @family External Cluster Validity Indices
//' @export
// [[Rcpp::export]]
double fowlkesMallows_fast(NumericVector pairCVIParams, bool normalize = false) {
  if (normalize) {
    return (pairCVIParams[0] - pairCVIParams[1]*pairCVIParams[2]/pairCVIParams[3]) /
      (sqrt(pairCVIParams[1]*pairCVIParams[2]) - pairCVIParams[1]*pairCVIParams[2]/pairCVIParams[3]);
  } else {
    return pairCVIParams[0] / sqrt(pairCVIParams[1]*pairCVIParams[2]);
  }
}

//' Phi Coefficient
//'
//' Calculates the Phi coefficient of Pearson (1900) to compare two cluster
//' assignment vectors (external cluster validation). It is a correlation value,
//' therefore being in the interval (-1,1], higher values indicating more similar
//' clusterings. The name varies in the literature. Instead of Phi, some sources
//' also call it Gamma. In Wu, Xiong and Chen (2009) its called "Hubert's Gamma
//' statistic I".
//'
//' @section References:
//'
//' Pearson, K. (1900). Mathematical contributions to the theory of evolution. vii.
//' on the correlation of characters not quantitatively measurable. \emph{Philosophical
//' Transactions of the Royal Society of London. Series A, Containing Papers of a
//' Mathematical or Physical Character}, 195, 1-405.
//'
//' Pearson, K. & Heron, D. (1913). On theories of association. \emph{Biometrika},
//' 9(1/2), 159-315.
//'
//' Wu, J., Xiong, H. & Chen, J. (2009). Adapting the right measures for k-means
//' clustering. In \emph{Proceedings of the 15th acm sigkdd international conference
//' on knowledge discovery and data mining} (pp. 877-886). ACM.
//'
//' @param pairCVIParams Output of \code{link{pairCVIParameters}} which has to be
//' called with the two cluster assignment vector to be compared.
//' @return Phi/Gamma coefficient as double from the range (-1,1].
//' @family External Cluster Validity Indices
//' @export
// [[Rcpp::export]]
double phi_fast(NumericVector pairCVIParams) {
  return (pairCVIParams[0]*pairCVIParams[3] - pairCVIParams[1]*pairCVIParams[2]) /
    sqrt(pairCVIParams[1] * pairCVIParams[2] * (pairCVIParams[3] - pairCVIParams[1]) *
      (pairCVIParams[3] - pairCVIParams[2]));
}

//' Purity Measure
//'
//' Calculates the purity measure (e.g. described by Wu, Xiong and Chen (2009))
//' to compare two cluster assignment vectors (external cluster validation). It
//' is a value in (0,1], higher values indicating more similarity. It finds the
//' most common ground truth class in each cluster and sums over these relative
//' frequencies.
//'
//' Be aware that this measure is asymmetric and can still be high if the classes
//' of the ground truth are split up into multiple (but pure) clusters. Wu, Xiong
//' and Chen (2009) propose to use the symmetric van Dongen measure
//' (\code{\link{vanDongen_fast}}) instead.
//'
//' @section References:
//'
//' Van Dongen, S. (2000). \emph{Performance criteria for graph clustering and markov
//' cluster experiments}. National Research Institute for Mathematics and Computer
//' Science. Amsterdam.
//'
//' Wu, J., Xiong, H. & Chen, J. (2009). Adapting the right measures for k-means
//' clustering. In \emph{Proceedings of the 15th acm sigkdd international conference
//' on knowledge discovery and data mining} (pp. 877-886). ACM.
//'
//' @param assignments Integer vector of cluster assignments containing only values
//' from 1 to k with k = number of clusters (code depends on this!).
//' @param groundTruth Integer vector of class (ground truth) assignments containing
//' only values from 1 to k with k = number of clusters.
//' @return The purity measure as double in (0,1].
//' @family External Cluster Validity Indices
//' @export
// [[Rcpp::export]]
double purity_fast(IntegerVector assignments, IntegerVector groundTruth) {
  IntegerMatrix countTable = IntegerMatrix(unique(assignments).size(), unique(groundTruth).size());
  for (int i = 0; i < assignments.size(); i++) {
    countTable(assignments[i] - 1, groundTruth[i] - 1)++;
  }
  double result = 0;
  for (int i = 0; i < countTable.nrow(); i++) {
    result += max(countTable.row(i));
  }
  result /= assignments.size(); // relative frequency
  return result;
}

//' Van Dongen Criterion
//'
//' Calculates the index of van Dongen (2000) to compare two cluster assignment
//' vectors (external cluster validation). It is a value in [0,2n), lower values
//' indicating more similarity (it matches each cluster of one assignment to the
//' most similar cluster in the other assignment and counts mismatches). Optionally,
//' the index can be normalized to [0,1] as proposed by Wu, Xiong and Chen (2009).
//' After normalization, we take 1 - normalizedValue so that higher values indicate
//' better clustering quality (as it is for indices like Rand, Fowlkes-Mallows).
//'
//' @section References:
//'
//' Van Dongen, S. (2000). \emph{Performance criteria for graph clustering and markov
//' cluster experiments}. National Research Institute for Mathematics and Computer
//' Science. Amsterdam.
//'
//' Wu, J., Xiong, H. & Chen, J. (2009). Adapting the right measures for k-means
//' clustering. In \emph{Proceedings of the 15th acm sigkdd international conference
//' on knowledge discovery and data mining} (pp. 877-886). ACM.
//'
//' @param assignments1 Integer vector of cluster assignments containing only values
//' from 1 to k with k = number of clusters (code depends on this!).
//' @param assignments2 Integer vector of cluster assignments containing only values
//' from 1 to k with k = number of clusters.
//' @param normalizeAndInvert Should the van Dongen criterion be normalized to
//' [0,1] and inverted such that high values indicate similar clusterings?
//' @return The van Dongen criterion as double (in [0,2n) without normalization
//' and [0,1] else).
//' @family External Cluster Validity Indices
//' @export
// [[Rcpp::export]]
double vanDongen_fast(IntegerVector assignments1, IntegerVector assignments2,
                      bool normalizeAndInvert = false) {
  IntegerVector table1 = IntegerVector(unique(assignments1).size());
  IntegerVector table2 = IntegerVector(unique(assignments2).size());
  IntegerMatrix table12 = IntegerMatrix(table1.size(), table2.size());
  for (int i = 0; i < assignments1.size(); i++) {
    table1[assignments1[i] - 1]++;
    table2[assignments2[i] - 1]++;
    table12(assignments1[i] - 1, assignments2[i] - 1)++;
  }
  int matches2to1 = 0, matches1to2 = 0;
  for (int i = 0; i < table12.nrow(); i++) {
    matches2to1 += max(table12.row(i));
  }
  for (int j = 0; j < table12.ncol(); j++) {
    matches1to2 += max(table12.column(j));
  }
  if (normalizeAndInvert) {
    // make sure that double division by multiplying with 2.0
    return 1 - (2.0 * assignments1.size() - matches1to2 - matches2to1) /
      (2.0 * assignments1.size() - max(table1) - max(table2));
  } else {
    return 2*assignments1.size() - matches1to2 - matches2to1;
  }
}

//' Variation of Information and Normalized Mutual Information
//'
//' Calculates the Variation of Information index introduced by of Meila (2003) to
//' compare two cluster assignment vectors (external cluster validation). It is a
//' value greater or equal 0, lower values indicating more similarity (it is based
//' on the entropy of the single assignments and the mutual information of the joint
//' distribution). Optionally, the index can be normalized to [0,1] as proposed by
//' Wu, Xiong and Chen (2009). After normalization, we take 1 - normalizedValue so
//' that higher values indicate better clustering quality (as it is for indices
//' like Rand, Fowlkes-Mallows); the result equals the Normalized Mutual Information
//' of Fred and Jain (2002).
//'
//' We use the base 2 logarithm for calculating entropy and mutual information.
//'
//' @section References:
//'
//' Fred, A. L. & Jain, A. K. (2002). Data clustering using evidence accumulation.
//' In \emph{Pattern recognition, 2002. proceedings. 16th international conference
//' on} (Vol. 4, pp. 276-280). IEEE.
//'
//' Meila, M. (2003). Comparing clusterings by the variation of information. In
//' B. Schölkopf & M. K. Warmuth (Eds.), \emph{Learning theory and kernel machines:
//' 16th annual conference on learning theory and 7th kernel workshop, colt/kernel
//' 2003, washington, dc, usa, august 24-27, 2003. proceedings} (pp. 173-187).
//' Springer Berlin Heidelberg.
//'
//' Wu, J., Xiong, H. & Chen, J. (2009). Adapting the right measures for k-means
//' clustering. In \emph{Proceedings of the 15th acm sigkdd international conference
//' on knowledge discovery and data mining} (pp. 877-886). ACM.
//'
//' @param assignments1 Integer vector of cluster assignments containing only values
//' from 1 to k with k = number of clusters (code depends on this!).
//' @param assignments2 Integer vector of cluster assignments containing only values
//' from 1 to k with k = number of clusters.
//' @param normalizeAndInvert Should the Variation of Information be normalized
//' to [0,1] and inverted such that high values indicate similar clusterings?
//' @return The Variation of Information as double (in [0, entropy1+entropy2] without
//' normalization and [0,1] else).
//' @export
// [[Rcpp::export]]
double VI_fast(IntegerVector assignments1, IntegerVector assignments2,
                                   bool normalizeAndInvert = false) {
  NumericVector table1 = NumericVector(unique(assignments1).size());
  NumericVector table2 = NumericVector(unique(assignments2).size());
  NumericMatrix table12 = NumericMatrix(table1.size(), table2.size());
  for (int i = 0; i < assignments1.size(); i++) {
    table1[assignments1[i] - 1]++;
    table2[assignments2[i] - 1]++;
    table12(assignments1[i] - 1, assignments2[i] - 1)++;
  }
  table1 = table1 / assignments1.size();
  table2 = table2 / assignments2.size();
  table12 = table12 / assignments1.size();
  double entropy1 = 0, entropy2 = 0, mutualInformation = 0;
  for (int i = 0; i < table1.size(); i++) {
    if (table1[i] != 0) {
      entropy1 -= table1[i] * log2(table1[i]);
    }
  }
  for (int j = 0; j < table2.size(); j++) {
    if (table2[j] != 0) {
      entropy2 -= table2[j] * log2(table2[j]);
    }
  }
  for (int i = 0; i < table12.nrow(); i++) {
    for (int j = 0; j < table12.ncol(); j++) {
      if (table12(i,j) != 0) {
        mutualInformation += table12(i,j) * log2(table12(i,j) / (table1[i] * table2[j]));
      }
    }
  }
  if (normalizeAndInvert) {
    return  2 * mutualInformation / (entropy1 + entropy2);
  } else {
    return entropy1 + entropy2 - 2 * mutualInformation;
  }
}
