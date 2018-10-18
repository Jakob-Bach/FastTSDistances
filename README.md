# Fast Dissimilarity Computations for Time Series

This R package contains fast (mostly C++) implementations of time series dissimilarities, simple aggregation functions for time series and cluster validity indices.
The package is not on `CRAN`, but can be installed manually.
Assuming you have installed the package `devtools`, you can simply run `devtools::install_github("Jakob-Bach/FastTSDistances")` to install `FastTSDistances` into your local R package repository.

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

- (external) Conditional entropy
- Entropy of an integer vector containing values from 1 to k (like cluster assignments)
- (external) Fowlkes-Mallows Index (optional: normalization)
- (internal) Generalized Davies-Bouldin Index based on intra- and inter-cluster dissimilarities; also inverted version which is higher for better validity
- (internal) Generalized Dunn Index based on intra- and inter-cluster dissimilarities
- (external) Phi Coefficient
- (external) Purity
- (external) Rand Index (optional: normalization)
- (external) van Dongen measure (optional: normalization)
- (external) Variation of Information (optional: normalization and inversion to Normalized Mutual Information)

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
