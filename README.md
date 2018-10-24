# Fast Dissimilarity Computations for Time Series

[![Build Status](https://travis-ci.com/Jakob-Bach/FastTSDistances.svg?branch=master)](https://travis-ci.com/Jakob-Bach/FastTSDistances)
[![codecov](https://codecov.io/gh/Jakob-Bach/FastTSDistances/branch/master/graph/badge.svg)](https://codecov.io/gh/Jakob-Bach/FastTSDistances)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This R package contains fast (mostly C++) implementations of time series dissimilarities, aggregation functions for time series and cluster validity indices.

## Overview

This package has been developed as part of a research project on understanding the effects of energy-data aggregation on clustering qualit.
For more information on the project, see the companion website: https://www.ipd.kit.edu/clustagg/.

## Setup
The FastTSDistances is not on `CRAN`, and has to be installed manually.
You can install the package directly from github via the [devtools package](https://github.com/r-lib/devtools).

```
devtools::install_github("Jakob-Bach/FastTSDistances")
```

After installation, call `help(package = "FastTSDistances")` for an overview of the package functions.

## Functionality

The package functionality falls into three categories: time-series dissimilarity, cluster validity, aggregation.

### Time-Series Dissimilarity

A dissimilarity function quantifies the dissimiliarity between two time series.
This package supports univariate time series (i.e., vectors) and multivariate time series (i.e., matrices with the columns representing different attributes).
It implements the following dissimilaritiy functions:

- Complexity-Invariant Distance (CID) for L2 metric and DTW
- Compression-/complexity-based dissimilarity (variant of CDM, using SAX representation and zip compression)
- Temporal correlation correction factor (CORT) for L2 metric and DTW
- Edit Distance on Real Sequences (EDR) (optional: Sakoe-Chiba window)
- Edit Distance with Real Penalty (ERP) (optional: Sakoe-Chiba window)
- Dynamic Time Warping (DTW) (optional: Sakoe-Chiba window)
- L1, L2, Lmax metric
- Permutation distribution dissimilarity
- Shaped-Based Distance (SBD)

### Cluster Validity

A cluster validty index quantifies the quality of a cluster assignment.
This package implements the following internal and external cluster validity indices:

* Entropy of an integer vector with cluster assignments from 1 to k
* External Validity
  * Conditional entropy (e.g. in [1])
  * Fowlkes-Mallows Index [2] (optional: with normalization [1])
  * Phi Coefficient [6]
  * Purity (e.g. in [1])
  * Rand Index [7, 8] (optional: normalization[1])
  * van Dongen measure [9] (optional: normalization [1])
  * Variation of Information [10] (optional: normalization [1] and inversion to Normalized Mutual Information [11])
* Internal Validity
  * Generalized Davies-Bouldin Index [3] based on intra- and inter-cluster dissimilarities; also in the inverted version which reports higher values for better validity
  * Generalized Dunn Index [4, 5] based on intra- and inter-cluster dissimilarities

### Aggregation

An aggregation function aggregates a univariate time series over fixed length intervals.
This package supports the following aggregation functions:

- piecewise kurtosis
- piecewise maximum
- piecewise mean
- piecewise median
- piecewise minimum
- piecewise skewness
- piecewise standard deviation
- Symbolic Aggregate Approximation (SAX)

### Further Dissimilarity-Related Functions

- (parallelized) dissimilarity matrix computation
- time series averaging
- z-scoring, min-max normalization


# References

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
