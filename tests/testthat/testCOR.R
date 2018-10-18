context("Correlation-based distance")

test_that("corDist_fast equivalent to TSclust::diss.COR (apart from some special cases)", {
  for (i in 1:50) {
    # not too short series, or precision difference between R implementation in
    # diss.COR and our own C++ implementation becomes too relevant for some
    # correlations (close to -1 and 1)
    n = floor(runif(n = 1, min = 10, max = 200))
    ts1 <- rnorm(n)
    ts2 <- rnorm(n)
    beta <- abs(rnorm(1))
    correlation <- cor(ts1, ts2)
    if (!is.na(correlation)) { # one series might be constant
      expect_equal(corDist_fast(ts1, ts2, beta),
                   TSclust::diss.COR(ts1 ,ts2, beta))
      expect_equal(corDist_fast(ts1, ts2, 0),
                   TSclust::diss.COR(ts1 ,ts2, NULL))
    }
  }
})

test_that("crossCorNormalized same as dtwclust::NCCc", {
  for (i in 1:50) {
    n = floor(runif(n = 1, min = 200, max = 1000))
    ts1 <- rnorm(n)
    ts2 <- rnorm(n)
    expect_equal(crossCorNormalized(ts1, ts2), dtwclust::NCCc(ts1 ,ts2))
  }
})

test_that("crossCorNormalized corner cases work", {
  for (i in 1:10) {
    n <- floor(runif(n = 1, min = 2, max = 200))
    ts1 <- rnorm(n)
    if (!all(ts1 == 0)) {
      expect_equal(crossCorNormalized(ts1, rep(0, n)), rep(0, 2*n - 1))
      expect_equal(crossCorNormalized(rep(0, n), ts1), rep(0, 2*n - 1))
      n2 <- floor(runif(n = 1, min = 2, max = 200))
      expect_equal(crossCorNormalized(rep(0, n), rep(0, n2)),
                   crossCorNormalized(rep(ts1[1], n), rep(ts1[1], n2)))
    }
  }
})

test_that("shapeBasedDistance (uni-variate) same as dtwclust::SBD", {
  for (i in 1:50) {
    n = floor(runif(n = 1, min = 2, max = 200))
    ts1 <- rnorm(n)
    ts2 <- rnorm(n)
    expect_equal(shapeBasedDistance(ts1,ts2), dtwclust::SBD(ts1,ts2)$dist)
  }
})

test_that("shapeBasedDistance (uni-variate) for constant series valid", {
  for (i in 1:10) {
    n = floor(runif(n = 1, min = 2, max = 200))
    ts1 <- rnorm(n)
    ts2 <- rep(rnorm(1), times = n)
    ccDist <- shapeBasedDistance(ts1,ts2)
    expect_true(ccDist >= 0 & ccDist <= 2)
  }
})

test_that("shapeBasedDistance (multi-variate) for identical attributes same as univariate", {
  for (i in 1:10) {
    n = floor(runif(n = 1, min = 2, max = 200))
    tsDim <- floor(runif(1, min = 2, max = 5))
    ts1 <- rnorm(n)
    ts2 <- rnorm(n)
    tsMult1 <- matrix(rep(ts1, times = tsDim), ncol = tsDim)
    tsMult2 <- matrix(rep(ts2, times = tsDim), ncol = tsDim)
    expect_equal(shapeBasedDistance(ts1, ts2), shapeBasedDistance(tsMult1, tsMult2))
  }
})

test_that("shapeBasedDistance (multi-variate) in expected range", {
  for (i in 1:50) {
    tsLength <- floor(runif(1, min = 2, max = 50))
    tsDim <- floor(runif(1, min = 2, max = 5))
    ts1 <- matrix(rnorm(n = tsLength * tsDim), nrow = tsLength)
    ts2 <- matrix(rnorm(n = tsLength * tsDim), nrow = tsLength)
    ccDist <- shapeBasedDistance(ts1,ts2)
    expect_true(ccDist >= 0 & ccDist <= 2)
  }
})

test_that("shapeBasedDistance (multi-variate) for constant series valid", {
  for (i in 1:10) {
    tsLength <- floor(runif(1, min = 2, max = 50))
    tsDim <- floor(runif(1, min = 2, max = 5))
    ts1 <- matrix(rnorm(n = tsLength * tsDim), nrow = tsLength)
    ts2 <- matrix(rnorm(n = tsLength * tsDim), nrow = tsLength)
    # Set random columns to a constant
    ts2[, sample(seq_len(tsDim), size = floor(runif(1, min = 1, max = tsDim + 1)))] <-
      rep(rnorm(1), times = tsLength)
    ccDist <- shapeBasedDistance(ts1,ts2)
    expect_true(ccDist >= 0 & ccDist <= 2)
  }
})
