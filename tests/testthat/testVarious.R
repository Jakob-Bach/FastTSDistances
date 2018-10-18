context("Various tests for remaining functions")

test_that("averageTimeSeries_fast is same as map-reduce implementation in R", {
  for (i in 1:50) {
    tsCount <- floor(runif(1, min = 1, max = 50))
    tsLength <- floor(runif(1, min = 1, max = 50))
    tsList <- lapply(seq_len(tsCount), function(i) rnorm(n = tsLength))
    tsWeights <- rnorm(tsCount)
    expect_equal(averageTimeSeries_fast(tsList, tsWeights),
                 Reduce("+", mapply("*", tsList, tsWeights, SIMPLIFY = FALSE)))
  }
})

test_that("averageTimeSeriesMult_fast is same as map-reduce implementation in R", {
  for (i in 1:50) {
    tsCount <- floor(runif(1, min = 1, max = 20))
    tsLength <- floor(runif(1, min = 1, max = 50))
    numAttributes <- floor(runif(1, min = 2, max = 10))
    tsList <- lapply(seq_len(tsCount), function(i) matrix(rnorm(tsLength * numAttributes),
                                                          ncol = numAttributes))
    tsWeights <- rnorm(tsCount)
    expect_equal(averageTimeSeriesMult_fast(tsList, tsWeights),
                 Reduce("+", mapply("*", tsList, tsWeights, SIMPLIFY = FALSE)))
  }
})

test_that("Elements of vector cross-distance is same as with outer()", {
  absDistance <- function(x,y) abs(x - y)
  for (i in 1:100) {
    ts1 <- rnorm(n = floor(runif(1, min = 1, max = 50)))
    ts2 <- rnorm(n = floor(runif(1, min = 1, max = 50)))
    distMat1 <- outer(ts1, ts2, FUN = absDistance)
    distMat2 <- vectorCrossDistMat(ts1, ts2)
    expect_equal(distMat2, distMat1)
  }
})

test_that("List of vector cross-distance is same as with proxy::dist", {
  for (i in 1:10) {
    tsLength <- floor(runif(1, min = 1, max = 50))
    tsList <- lapply(1:50, function(i) rnorm(n = tsLength))
    expect_equal(sum(abs(proxy::as.matrix(proxy::dist(tsList, method = "L2")) -
                           tsCrossDistMat(tsList, distMethod = "l2Dist_fast", trace = FALSE)))
                 < 1e-10, TRUE)
  }
})

test_that("minMaxNormalize is in [0,1]", {
  for (i in 1:50) {
    ts <- rnorm(n = floor(runif(1, min = 2, max = 100)))
    maxIdx <- which.max(ts)
    minIdx <- which.min(ts)
    tsNorm <- minMaxNormalize(ts)
    expect_true(all(tsNorm >= 0 & tsNorm <= 1))
    expect_equal(tsNorm[maxIdx], 1)
    expect_equal(tsNorm[minIdx], 0)
  }
})

test_that("minMaxNormalize multivariate is same as minMaxNormalize per column", {
  for (i in 1:50) {
    numAttributes <- floor(runif(1, min = 2, max = 10))
    n <- floor(runif(1, min = 2, max = 100))
    ts <- matrix(rnorm(n * numAttributes), ncol = numAttributes)
    tsNorm <- minMaxNormalize(ts)
    tsNormExpected <- apply(ts, 2, minMaxNormalize)
    expect_equal(tsNorm, tsNormExpected)
  }
})

test_that("znormalize is same as proxy::base", {
  for (i in 1:50) {
    ts <- rnorm(n = floor(runif(1, min = 2, max = 100)))
    if (sd(ts) != 0) {
      expect_equal(znormalize(ts), as.numeric(base::scale(ts)))
    }
  }
})

test_that("znormalize multivariate is same as proxy::base per column", {
  for (i in 1:50) {
    numAttributes <- floor(runif(1, min = 2, max = 10))
    n <- floor(runif(1, min = 2, max = 100))
    ts <- matrix(rnorm(n * numAttributes), ncol = numAttributes)
    if (all(apply(ts, 2, sd) != 0)) {
      tsNorm <- znormalize(ts)
      tsNormExpected <- apply(ts, 2, base::scale)
      expect_equal(tsNorm, tsNormExpected)
    }
  }
})

test_that("znormalize can handle constant ts", {
  for (i in 1:10) {
    ts <- rep(rnorm(1), times = floor(runif(1, min = 1, max = 100)))
    expect_equal(znormalize(ts), rep(0, times = length(ts)))
  }
})
