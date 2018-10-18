context("L_p norm/metric")

test_that("l1Dist_fast is same as proxy's L1 distance", {
  for (i in 1:50) {
    tsLength <- floor(runif(1, min = 1, max = 50))
    ts1 <- rnorm(n = tsLength)
    ts2 <- rnorm(n = tsLength)
    expect_equal(l1Dist_fast(ts1, ts2),
                 as.numeric(proxy::dist(list(ts1), list(ts2), method = "L1", pairwise = TRUE)))
  }
})

test_that("l1DistMult_fast is same as with manual computation", {
  for (i in 1:50) {
    tsLength <- floor(runif(1, min = 1, max = 50))
    tsDim <- floor(runif(1, min = 2, max = 5))
    ts1 <- matrix(rnorm(n = tsLength * tsDim), nrow = tsLength)
    ts2 <- matrix(rnorm(n = tsLength * tsDim), nrow = tsLength)
    expect_equal(l1DistMult_fast(ts1, ts2),
                 sum(abs(ts1 - ts2)))
  }
})

test_that("l2Dist is same as proxy's L2 distance", {
  for (i in 1:50) {
    tsLength <- floor(runif(1, min = 1, max = 50))
    ts1 <- rnorm(n = tsLength)
    ts2 <- rnorm(n = tsLength)
    expect_equal(l2Dist(ts1, ts2),
                 as.numeric(proxy::dist(list(ts1), list(ts2), method = "L2", pairwise = TRUE)))
  }
})

test_that("l2Dist_fast is same as proxy's L2 distance", {
  for (i in 1:50) {
    tsLength <- floor(runif(1, min = 1, max = 50))
    ts1 <- rnorm(n = tsLength)
    ts2 <- rnorm(n = tsLength)
    expect_equal(l2Dist_fast(ts1, ts2),
                 as.numeric(proxy::dist(list(ts1), list(ts2), method = "L2", pairwise = TRUE)))
  }
})

test_that("l2DistMult_fast is same as with manual computation", {
  for (i in 1:50) {
    tsLength <- floor(runif(1, min = 1, max = 50))
    tsDim <- floor(runif(1, min = 2, max = 5))
    ts1 <- matrix(rnorm(n = tsLength * tsDim), nrow = tsLength)
    ts2 <- matrix(rnorm(n = tsLength * tsDim), nrow = tsLength)
    expect_equal(l2DistMult_fast(ts1, ts2),
                 sqrt(sum((ts1 - ts2)^2)))
  }
})

test_that("l2Dist_fast with cid is same as TSclust::diss.CID", {
  for (i in 1:50) {
    tsLength <- floor(runif(1, min = 2, max = 50))
    ts1 <- rnorm(n = tsLength)
    ts2 <- rnorm(n = tsLength)
    expect_equal(l2Dist_fast(ts1, ts2, cid = TRUE),
                 as.numeric(TSclust::diss.CID(ts1,ts2)))
  }
})

test_that("l2Dist_fast with cid is l2Dist_fast for two constant series", {
  for (i in 1:10) {
    tsLength <- floor(runif(1, min = 2, max = 50))
    ts1 <- rep(rnorm(1), times = tsLength)
    ts2 <- rep(rnorm(1), times = tsLength)
    expect_equal(l2Dist_fast(ts1, ts2, cid = TRUE),
                 l2Dist_fast(ts1, ts2, cid = FALSE))
  }
})

test_that("l2DistMult_fast with cid is same as with manual computation", {
  for (i in 1:50) {
    tsLength <- floor(runif(1, min = 2, max = 50))
    tsDim <- floor(runif(1, min = 2, max = 5))
    ts1 <- matrix(rnorm(n = tsLength * tsDim), nrow = tsLength)
    ts2 <- matrix(rnorm(n = tsLength * tsDim), nrow = tsLength)
    complexityX <- sqrt(sum(apply(ts1, 2, function(x) sum(diff(x)^2))))
    complexityY <- sqrt(sum(apply(ts2, 2, function(x) sum(diff(x)^2))))
    expect_equal(l2DistMult_fast(ts1, ts2, cid = TRUE),
                 l2DistMult_fast(ts1, ts2) * max(complexityX, complexityY) /
                   min(complexityX, complexityY))
  }
})

test_that("l2DistMult_fast with cid is l2DistMult_fast for two series with constant attributes", {
  for (i in 1:10) {
    tsLength <- floor(runif(1, min = 2, max = 50))
    tsDim <- floor(runif(1, min = 2, max = 5))
    ts1 <- sapply(seq_along(tsDim), function(i) rep(rnorm(1), times = tsLength))
    ts2 <- sapply(seq_along(tsDim), function(i) rep(rnorm(1), times = tsLength))
    expect_equal(l2Dist_fast(ts1, ts2, cid = TRUE),
                 l2Dist_fast(ts1, ts2, cid = FALSE))
  }
})

test_that("l2Dist_fast with cort is same as TSclust::diss.CORT", {
  for (i in 1:50) {
    tsLength <- floor(runif(1, min = 2, max = 50))
    ts1 <- rnorm(n = tsLength)
    ts2 <- rnorm(n = tsLength)
    k <- abs(runif(1))
    expect_equal(l2Dist_fast(ts1, ts2, cortK = k),
                 as.numeric(TSclust::diss.CORT(ts1,ts2, k = k, deltamethod = "Euclid")))
  }
})

test_that("l2Dist_fast with cort also works for constant series", {
  for (i in 1:20) {
    tsLength <- floor(runif(1, min = 2, max = 50))
    ts1 <- rep(rnorm(1), times = tsLength)
    k <- abs(runif(1))
    if (i %% 2 == 1) {
      ts2 <- rep(rnorm(1), times = tsLength)
      expect_equal(l2Dist_fast(ts1, ts2, cortK = k), 2 / (1 + exp(k)) * l2Dist_fast(ts1, ts2))
    } else {
      ts2 <- rnorm(n = tsLength)
      expect_equal(l2Dist_fast(ts1, ts2, cortK = k), l2Dist_fast(ts1, ts2))
    }
  }
})

test_that("l2DistMult_fast with cort is same as TSclust::diss.CORT for 1D matrix", {
  # implementation in TSclust is only univariate
  for (i in 1:50) {
    tsLength <- floor(runif(1, min = 2, max = 50))
    ts1 <- rnorm(n = tsLength)
    ts2 <- rnorm(n = tsLength)
    k <- abs(runif(1))
    expect_equal(l2DistMult_fast(as.matrix(ts1), as.matrix(ts2), cortK = k),
                 as.numeric(TSclust::diss.CORT(ts1,ts2, k = k, deltamethod = "Euclid")))
  }
})

test_that("cortFactorMult_fast for identical attributes is same as cortFactor_fast", {
  for (i in 1:10) {
    n = floor(runif(n = 1, min = 2, max = 200))
    tsDim <- floor(runif(1, min = 2, max = 5))
    ts1 <- rnorm(n)
    ts2 <- rnorm(n)
    tsMult1 <- matrix(rep(ts1, times = tsDim), ncol = tsDim)
    tsMult2 <- matrix(rep(ts2, times = tsDim), ncol = tsDim)
    expect_equal(cortFactor_fast(ts1, ts2), cortFactorMult_fast(tsMult1, tsMult2))
  }
})

test_that("cortFactorMult_fast is in [0,2]", {
  for (i in 1:50) {
    tsLength <- floor(runif(1, min = 2, max = 50))
    tsDim <- floor(runif(1, min = 2, max = 5))
    ts1 <- matrix(rnorm(n = tsLength * tsDim), nrow = tsLength)
    if (i %% 10 == 1) {
      # Set random columns to a constant
      ts1[, sample(seq_len(tsDim), size = floor(runif(1, min = 1, max = tsDim + 1)))] <-
        rep(rnorm(1), times = tsLength)
    }
    ts2 <- matrix(rnorm(n = tsLength * tsDim), nrow = tsLength)
    k <- abs(runif(1))
    cortFactor <- cortFactorMult_fast(ts1, ts2, k = k)
    expect_true(cortFactor >= 0 && cortFactor <= 2)
  }
})

test_that("lmaxDist_fast is same as proxy's Lmax distance", {
  for (i in 1:50) {
    tsLength <- floor(runif(1, min = 1, max = 50))
    ts1 <- rnorm(n = tsLength)
    ts2 <- rnorm(n = tsLength)
    expect_equal(lmaxDist_fast(ts1, ts2),
                 as.numeric(proxy::dist(list(ts1), list(ts2), method = "Chebyshev", pairwise = TRUE)))
  }
})

test_that("lmaxDistMult_fast is same as with manual computation", {
  for (i in 1:50) {
    tsLength <- floor(runif(1, min = 1, max = 50))
    tsDim <- floor(runif(1, min = 2, max = 5))
    ts1 <- matrix(rnorm(n = tsLength * tsDim), nrow = tsLength)
    ts2 <- matrix(rnorm(n = tsLength * tsDim), nrow = tsLength)
    expect_equal(lmaxDistMult_fast(ts1, ts2),
                 max(abs(ts1 - ts2)))
  }
})

test_that("l2Norm_fast acually is l2 norm", {
  for (i in 1:50) {
    tsLength <- floor(runif(1, min = 1, max = 50))
    ts <- rnorm(n = tsLength)
    expect_equal(l2Norm_fast(ts), sqrt(sum(ts*ts)))
  }
})
