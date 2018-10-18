context("Edit distances")

test_that("EDRDist_fast equivalent to TSDist::EDRDistance", {
  for (i in 1:50) {
    ts1 <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    ts2 <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    epsilon <- abs(rnorm(1))
    expect_equal(EDRDist_fast(ts1, ts2, epsilon), TSdist::EDRDistance(ts1, ts2, epsilon))
  }
})

test_that("EDRDistMult_fast equivalent to TSDist::EDRDistance for 1D matrix", {
  # implementation in TSdist is only univariate
  for (i in 1:50) {
    ts1 <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    ts2 <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    epsilon <- abs(rnorm(1))
    expect_equal(EDRDistMult_fast(as.matrix(ts1), as.matrix(ts2), epsilon),
                 TSdist::EDRDistance(ts1, ts2, epsilon))
  }
})

test_that("EDRDistSakoeChiba_fast equivalent to TSDist::EDRDistance", {
  for (i in 1:50) {
    ts1 <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    ts2 <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    epsilon <- abs(rnorm(1))
    # No solution if window smaller than length difference
    minWindowSize <- max(abs(length(ts1) - length(ts2)), 1)
    windowSize <- floor(runif(n = 1, min = minWindowSize, max = max(length(ts1), length(ts2))))
    expect_equal(EDRDistSakoeChiba_fast(ts1, ts2, epsilon, windowSize),
                 max(TSdist::EDRDistance(ts1, ts2, epsilon, windowSize),
                     abs(length(ts1) - length(ts2))))
    # There seems to be a bug in TSdist::ERPDistance(): if the window size is exactly the
    # length difference + all pairwise distances are under the threshold epsilon, then 1
    # is returned instead of the length difference which would be correct. As a consequence,
    # we make sure that the expected value is at least the length difference by using max()
  }
})

test_that("EDRDistSakoeChibaMult_fast equivalent to TSDist::EDRDistance for 1D matrix", {
  # implementation in TSdist is only univariate
  for (i in 1:50) {
    ts1 <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    ts2 <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    epsilon <- abs(rnorm(1))
    # No solution if window smaller than length difference
    minWindowSize <- max(abs(length(ts1) - length(ts2)), 1)
    windowSize <- floor(runif(n = 1, min = minWindowSize, max = max(length(ts1), length(ts2))))
    expect_equal(EDRDistSakoeChibaMult_fast(as.matrix(ts1), as.matrix(ts2), epsilon, windowSize),
                 max(TSdist::EDRDistance(ts1, ts2, epsilon, windowSize),
                     abs(length(ts1) - length(ts2)))) # see uni-variate routine
  }
})

test_that("ERPDist equivalent to TSDist::ERPDistance", {
  for (i in 1:50) {
    ts1 <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    ts2 <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    gapPenalty <- rnorm(1)
    expect_equal(ERPDist(ts1, ts2, gapPenalty), TSdist::ERPDistance(ts1, ts2, gapPenalty))
  }
})

test_that("ERPDist_fast equivalent to TSDist::ERPDistance", {
  for (i in 1:50) {
    ts1 <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    ts2 <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    gapPenalty <- rnorm(1)
    expect_equal(ERPDist_fast(ts1, ts2, gapPenalty), TSdist::ERPDistance(ts1, ts2, gapPenalty))
  }
})

test_that("ERPDistMult_fast equivalent to TSDist::ERPDistance for 1D matrix", {
  # implementation in TSdist is only univariate
  for (i in 1:50) {
    ts1 <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    ts2 <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    gapPenalty <- rnorm(1)
    expect_equal(ERPDistMult_fast(as.matrix(ts1), as.matrix(ts2), gapPenalty),
                 TSdist::ERPDistance(ts1, ts2, gapPenalty))
  }
})

test_that("ERPDistSakoeChiba equivalent to TSDist::ERPDistance", {
  for (i in 1:50) {
    ts1 <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    ts2 <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    gapPenalty <- rnorm(1)
    # No solution if window smaller than length difference
    minWindowSize <- max(abs(length(ts1) - length(ts2)), 1)
    windowSize <- floor(runif(n = 1, min = minWindowSize, max = max(length(ts1), length(ts2))))
    expect_equal(ERPDistSakoeChiba(ts1, ts2, gapPenalty, windowSize),
                 TSdist::ERPDistance(ts1, ts2, gapPenalty, windowSize))
  }
})

test_that("ERPDistSakoeChiba_fast equivalent to TSDist::ERPDistance", {
  for (i in 1:50) {
    ts1 <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    ts2 <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    gapPenalty <- rnorm(1)
    # No solution if window smaller than length difference
    minWindowSize <- max(abs(length(ts1) - length(ts2)), 1)
    windowSize <- floor(runif(n = 1, min = minWindowSize, max = max(length(ts1), length(ts2))))
    expect_equal(ERPDistSakoeChiba_fast(ts1, ts2, gapPenalty, windowSize),
                 TSdist::ERPDistance(ts1, ts2, gapPenalty, windowSize))
  }
})

test_that("ERPDistSakoeChibaMult_fast equivalent to TSDist::ERPDistance", {
  # implementation in TSdist is only univariate
  for (i in 1:50) {
    ts1 <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    ts2 <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    gapPenalty <- rnorm(1)
    # No solution if window smaller than length difference
    minWindowSize <- max(abs(length(ts1) - length(ts2)), 1)
    windowSize <- floor(runif(n = 1, min = minWindowSize, max = max(length(ts1), length(ts2))))
    expect_equal(ERPDistSakoeChibaMult_fast(as.matrix(ts1), as.matrix(ts2), gapPenalty, windowSize),
                 TSdist::ERPDistance(ts1, ts2, gapPenalty, windowSize))
  }
})
