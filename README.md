# Fast Dissimilarity Computations for Time Series

[![Build Status](https://travis-ci.com/Jakob-Bach/FastTSDistances.svg?branch=master)](https://travis-ci.com/Jakob-Bach/FastTSDistances)
[![codecov](https://codecov.io/gh/Jakob-Bach/FastTSDistances/branch/master/graph/badge.svg)](https://codecov.io/gh/Jakob-Bach/FastTSDistances)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This R package contains fast (mostly C++) implementations of time series dissimilarities, simple aggregation functions for time series and cluster validity indices.
The package is not on `CRAN`, but can be installed manually.
Assuming you have installed the package `devtools`, you can simply run `devtools::install_github("Jakob-Bach/FastTSDistances")` to install `FastTSDistances` into your local R package repository.
After installation, you can call `help(package = "FastTSDistances")` to get an overview of the package's functions.

## Aggregation Functions

All aggregation functions aggregate intervals of a fixed length from univariate time series (=vectors).

- piecewise kurtosis
- piecewise maximum
- piecewise mean [1, 2]
- piecewise median
- piecewise minimum
- piecewise skewness
- piecewise standard deviation
- Symbolic Aggregate Approximation [3] (SAX)

[1] Keogh, E. J. & Pazzani, M. J. (2000). Scaling up dynamic time warping for datamining applications.

[2] Keogh, E., Chakrabarti, K., Pazzani, M. & Mehrotra, S. (2001). Dimensionality reduction for fast similarity search in large time series databases.

[3] Lin, J., Keogh, E., Lonardi, S. & Chiu, B. (2003). A symbolic representation of time series, with implications for streaming algorithms.

## Cluster Validity Indices

We provide internal as well as external cluster validity indices.

- (external) Conditional entropy (e.g. in [1])
- Entropy of an integer vector containing values from 1 to k (like cluster assignments)
- (external) Fowlkes-Mallows Index [2] (optional: normalization [1])
- (internal) Generalized Davies-Bouldin Index [3] based on intra- and inter-cluster dissimilarities; also inverted version which is higher for better validity
- (internal) Generalized Dunn Index [4, 5] based on intra- and inter-cluster dissimilarities
- (external) Phi Coefficient [6]
- (external) Purity (e.g. in [1])
- (external) Rand Index [7, 8] (optional: normalization[1])
- (external) van Dongen measure [9] (optional: normalization [1])
- (external) Variation of Information [10] (optional: normalization [1] and inversion to Normalized Mutual Information [11])

[1] Wu, J., Xiong, H. & Chen, J. (2009). Adapting the right measures for k-means clustering.

[2] Fowlkes, E. B. & Mallows, C. L. (1983). A method for comparing two hierarchical clusterings.

[3] Davies, D. L. & Bouldin, D. W. (1979). A cluster separation measure.

[4] Dunn, J. C. (1973). A fuzzy relative of the isodata process and its use in detecting compact well-separated clusters.

[5] Bezdek, J. C. & Pal, N. R. (1998). Some new indexes of cluster validity.

[6] Pearson, K. (1900). Mathematical contributions to the theory of evolution. vii. on the correlation of characters not quantitatively measurable.

[7] Rand, W. M. (1971). Objective criteria for the evaluation of clustering methods.

[8] Hubert, L. & Arabie, P. (1985). Comparing partitions.

[9] Van Dongen, S. (2000). Performance criteria for graph clustering and markov cluster experiments.

[10] Meila, M. (2003). Comparing clusterings by the variation of information.

[11] Fred, A. L. & Jain, A. K. (2002). Data clustering using evidence accumulation.

## Dissimilarity Functions

All dissimilarities support univariate time series (=vectors) as well as multivariate time series (matrices with the columns representing different attributes).

- Complexity-Invariant Distance (CID) [1] for L2 metric and DTW
- Compression-/complexity-based dissimilarity [2, 3] (variant of CDM, using SAX representation and zip compression)
- Correlaion-based distance [4]
- Temporal correlation correction factor (CORT) [5] for L2 metric and DTW
- Edit Distance on Real Sequences (EDR) [6] (optional: Sakoe-Chiba window [7])
- Edit Distance with Real Penalty (ERP) [8] (optional: Sakoe-Chiba window [7])
- Dynamic Time Warping (DTW) [9] (optional: Sakoe-Chiba window [7])
- L1, L2, Lmax metric
- Permutation distribution dissimilarity [10]
- Shaped-Based Distance (SBD) [11]

[1] Batista, G. E., Keogh, E. J., Tataw, O. M. & De Souza, V. M. (2014). Cid: An efficient complexity-invariant distance for time series.

[2] Li, M., Badger, J. H., Chen, X., Kwong, S., Kearney, P. & Zhang, H. (2001). An information-based sequence distance and its application to whole mitochondrial genome phylogeny.

[3] Keogh, E., Lonardi, S., Ratanamahatana, C. A., Wei, L., Lee, S.-H. & Handley, J. (2007). Compression-based data mining of sequential data.

[4] Golay, X., Kollias, S., Stoll, G., Meier, D., Valavanis, A. & Boesiger, P. (1998). A new correlation-based fuzzy logic clustering algorithm for fmri.

[5] Chouakria, A. D. & Nagabhushan, P. N. (2007). Adaptive dissimilarity index for measuring time series proximity.

[6] Chen, L., Özsu, M. T. & Oria, V. (2005). Robust and fast similarity search for moving object trajectories.

[7] Sakoe, H., & Chiba, S. (1978). Dynamic programming algorithm optimization for spoken word recognition.

[8] Chen, L., & Ng, R. (2004, August). On the marriage of lp-norms and edit distance.

[9] Berndt, D. J. & Clifford, J. (1994). Using dynamic time warping to find patterns in time series.

[10] Brandmaier, A. M. (2011). Permutation distribution clustering and structural equation model trees.

[11] Paparrizos, J. & Gravano, L. (2015). K-shape: Efficient and accurate clustering of time series.

## General Dissimilarity-Related Functions

- (parallelized) dissimilarity matrix computation
- time series averaging
- z-scoring, min-max normalization
