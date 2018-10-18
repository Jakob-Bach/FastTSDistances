context("Aggregation")

test_that("PAA_fast equivalent to TSclust::PAA", {
  for (i in 1:10) {
    ts <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    for (windowCount in c(1,length(ts), sample(x = 1:length(ts), size = 8, replace = TRUE))) {
      expect_equal(PAA_fast(ts, windowCount),
                   TSclust::PAA(ts, windowCount))
    }
  }
})

test_that("SAX_fast/SAXLimitsOriginal equivalent to TSclust::convert.to.SAX.symbol", {
  for (i in 1:10) {
    ts <- rnorm(n = floor(runif(n = 1, min = 2, max = 200)))
    for (symbolCount in 1:10) {
      expect_equal(SAX_fast(ts, SAXLimitsOriginal(symbolCount)),
                   TSclust::convert.to.SAX.symbol(ts, symbolCount))
    }
  }
})

segmentToList <- function(x, segmentCount) {
  segmentLength <- length(x) / segmentCount
  return(lapply(1:segmentCount, function(i) {
    x[floor((i - 1) * segmentLength + 1):floor(i * segmentLength)]
  }))
}

test_that("PMaxAA_fast equivalent to piecewise max", {
  for (i in 1:100) {
    segmentLength <- floor(runif(n = 1, min = 1, max = 20))
    segmentCount <- floor(runif(n = 1, min = 1, max = 10))
    ts <- rnorm(n = segmentLength * segmentCount)
    expect_equal(PMaxAA_fast(ts, segmentCount), sapply(segmentToList(ts, segmentCount), max))
  }
})

test_that("PMedAA_fast equivalent to piecewise median", {
  for (i in 1:100) {
    segmentLength <- floor(runif(n = 1, min = 1, max = 20))
    segmentCount <- floor(runif(n = 1, min = 1, max = 10))
    ts <- rnorm(n = segmentLength * segmentCount)
    expect_equal(PMedAA_fast(ts, segmentCount), sapply(segmentToList(ts, segmentCount), median))
  }
})

test_that("PMinAA_fast equivalent to piecewise min", {
  for (i in 1:100) {
    segmentLength <- floor(runif(n = 1, min = 1, max = 20))
    segmentCount <- floor(runif(n = 1, min = 1, max = 10))
    ts <- rnorm(n = segmentLength * segmentCount)
    expect_equal(PMinAA_fast(ts, segmentCount), sapply(segmentToList(ts, segmentCount), min))
  }
})

test_that("PSDAA_fast equivalent to piecewise sd", {
  for (i in 1:50) {
    segmentLength <- floor(runif(n = 1, min = 1, max = 20))
    segmentCount <- floor(runif(n = 1, min = 1, max = 10))
    ts <- rnorm(n = segmentLength * segmentCount)
    expect_equal(PSDAA_fast(ts, segmentCount),
                 sapply(segmentToList(ts, segmentCount), function(segment) {
                   Weighted.Desc.Stat::w.sd(segment, rep(1, length(segment)))
                 }))
  }
})

test_that("PSDAA_fast (sample) equivalent to piecewise sd", {
  for (i in 1:50) {
    segmentLength <- floor(runif(n = 1, min = 1, max = 20))
    segmentCount <- floor(runif(n = 1, min = 1, max = 10))
    ts <- rnorm(n = segmentLength * segmentCount)
    expect_equal(PSDAA_fast(ts, segmentCount, sample = TRUE),
                 sapply(segmentToList(ts, segmentCount), sd))
  }
})

test_that("PSkewAA_fast equivalent to piecewise skewness", {
  for (i in 1:100) {
    # min segment length of two, because otherwise division 0/0 which might not
    # always be NaN in C++ (is undefined, might even differ on same machine if
    # both i386 and x64 tests are run during devtools::check)
    segmentLength <- floor(runif(n = 1, min = 2, max = 20))
    segmentCount <- floor(runif(n = 1, min = 1, max = 10))
    ts <- rnorm(n = segmentLength * segmentCount)
    expect_equal(PSkewAA_fast(ts, segmentCount),
                 sapply(segmentToList(ts, segmentCount), function(segment) {
                   Weighted.Desc.Stat::w.skewness(segment, rep(1, length(segment)))
                 }))
  }
})

test_that("PKurtAA_fast equivalent to piecewise kurtosis", {
  for (i in 1:100) {
    # min segment length of two, because otherwise division 0/0 which might not
    # always be NaN in C++ (is undefined, might even differ on same machine if
    # both i386 and x64 tests are run during devtools::check)
    segmentLength <- floor(runif(n = 1, min = 2, max = 20))
    segmentCount <- floor(runif(n = 1, min = 1, max = 10))
    ts <- rnorm(n = segmentLength * segmentCount)
    expect_equal(PKurtAA_fast(ts, segmentCount),
                 sapply(segmentToList(ts, segmentCount), function(segment) {
                   # w.kurtosis computes excess kurtosis, which is kurtosis - 3
                   Weighted.Desc.Stat::w.kurtosis(segment, rep(1, length(segment))) + 3
                 }))
  }
})
