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
- piecewise mean
- piecewise median
- piecewise minimum
- piecewise skewness
- piecewise standard deviation
- Symbolic Aggregate Approximation (SAX)

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

- Complexity-Invariant Distance (CID) for L2 metric and DTW
- Compression-/complexity-based dissimilarity (variant of CDM, using SAX representation and zip compression)
- Temporal correlation correction factor (CORT) for L2 metric and DTW
- Edit Distance on Real Sequences (EDR) (optional: Sakoe-Chiba window)
- Edit Distance with Real Penalty (ERP) (optional: Sakoe-Chiba window)
- Dynamic Time Warping (DTW) (optional: Sakoe-Chiba window)
- L1, L2, Lmax metric
- Permutation distribution dissimilarity
- Shaped-Based Distance (SBD)

## General Dissimilarity-Related Functions

- (parallelized) dissimilarity matrix computation
- time series averaging
- z-scoring, min-max normalization
