context("DTW")

test_that("DTWDist_fast equivalent to dtwclust::dtw_basic", {
  for (i in 1:50) {
    ts1 <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    ts2 <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    expect_equal(DTWDist_fast(ts1, ts2),
                 dtw::dtw(ts1, ts2, step.pattern = dtw::symmetric1,
                          distance.only = TRUE)$distance)
  }
})

test_that("DTWDistMult_fast equivalent to dtwclust::dtw_basic", {
  for (i in 1:50) {
    tsLength <- floor(runif(1, min = 1, max = 50))
    tsDim <- floor(runif(1, min = 2, max = 5))
    ts1 <- matrix(rnorm(n = tsLength * tsDim), nrow = tsLength)
    ts2 <- matrix(rnorm(n = tsLength * tsDim), nrow = tsLength)
    expect_equal(DTWDistMult_fast(ts1, ts2),
                 dtw::dtw(ts1, ts2, step.pattern = dtw::symmetric1,
                          distance.only = TRUE)$distance)
  }
})

# test_that("DTWDist_fast with cort equivalent to TSclust::diss:CORT", {
#   # no test here as dtw::dtw() which used in diss.CORT() has a different step
#   # pattern compared to our DTW and we cannot configure it in the method
#   # call of diss.CORT()
# })

test_that("DTWDistSakoeChiba_fast equivalent to dtwclust::dtw_basic", {
  for (i in 1:50) {
    ts1 <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    ts2 <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    # No solution if window smaller than length difference
    minWindowSize <- max(abs(length(ts1) - length(ts2)), 1)
    windowSize <- floor(runif(n = 1, min = minWindowSize, max = max(length(ts1), length(ts2))))
    expect_equal(DTWDistSakoeChiba_fast(ts1, ts2, windowSize),
                 dtw::dtw(ts1, ts2, step.pattern = dtw::symmetric1, distance.only = TRUE,
                          window.type = "sakoechiba", window.size = windowSize)$distance)
  }
})

test_that("DTWDistSakoeChibaMult_fast equivalent to dtwclust::dtw_basic", {
  for (i in 1:50) {
    tsLength <- floor(runif(1, min = 1, max = 50))
    tsDim <- floor(runif(1, min = 2, max = 5))
    ts1 <- matrix(rnorm(n = tsLength * tsDim), nrow = tsLength)
    ts2 <- matrix(rnorm(n = tsLength * tsDim), nrow = tsLength)
    # No solution if window smaller than length difference
    minWindowSize <- max(abs(length(ts1) - length(ts2)), 1)
    windowSize <- floor(runif(n = 1, min = minWindowSize, max = max(length(ts1), length(ts2))))
    expect_equal(DTWDistSakoeChibaMult_fast(ts1, ts2, windowSize),
                 dtw::dtw(ts1, ts2, step.pattern = dtw::symmetric1, distance.only = TRUE,
                          window.type = "sakoechiba", window.size = windowSize)$distance)
  }
})
