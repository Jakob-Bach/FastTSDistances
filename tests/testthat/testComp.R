context("Complexity-/compression-based dissimilarities")

test_that("compDist(x,x) (uni-variate) is 0", {
  for (i in 1:10) {
    ts <- rnorm(n = floor(runif(1, min = 1, max = 200)))
    expect_equal(compDist(ts,ts), 0)
  }
})

test_that("compDist(x,x) (multi-variate) is 0", {
  for (i in 1:10) {
    tsLength <- floor(runif(1, min = 1, max = 50))
    tsDim <- floor(runif(1, min = 2, max = 5))
    ts <- matrix(rnorm(n = tsLength * tsDim), nrow = tsLength)
    expect_equal(compDist(ts,ts), 0)
  }
})

test_that("compDist(x,y) (uni-variate) is in [0,1] and symmetric", {
  for (i in 1:50) {
    saxLength <- floor(runif(1, min = 1, max = 20))
    tsLength <- floor(runif(1, min = 1, max = 50))
    tsDim <- floor(runif(1, min = 2, max = 5))
    ts1 <- matrix(rnorm(n = tsLength * tsDim), nrow = tsLength)
    ts2 <- matrix(rnorm(n = tsLength * tsDim), nrow = tsLength)
    distance12 <- compDist(ts1, ts2, symbolCount = saxLength)
    distance21 <- compDist(ts2, ts1, symbolCount = saxLength)
    expect_true(distance12 >= 0 && distance12 <= 1)
    expect_true(distance21 >= 0 && distance21 <= 1)
    expect_equal(distance12, distance21)
  }
})

test_that("compDist(x,y) (multi-variate) is in [0,1] and symmetric", {
  for (i in 1:50) {
    tsLength <- floor(runif(1, min = 1, max = 50))
    tsDim <- floor(runif(1, min = 2, max = 5))
    saxLength <- floor(runif(1, min = 1, max = 20))
    ts1 <- matrix(rnorm(tsLength * tsDim), ncol = tsDim)
    ts2 <- matrix(rnorm(tsLength * tsDim), ncol = tsDim)
    distance12 <- compDist(ts1, ts2, symbolCount = saxLength)
    distance21 <- compDist(ts2, ts1, symbolCount = saxLength)
    expect_true(distance12 >= 0 && distance12 <= 1)
    expect_true(distance21 >= 0 && distance21 <= 1)
    expect_equal(distance12, distance21)
  }
})

test_that("compDistList same as nested loop calls", {
  for (i in 1:16) {
    tsLength <- floor(runif(1, min = 1, max = 50))
    symbolCount <- floor(runif(1, min = 1, max = 10))
    if (i %% 4 == 1) {# uni-variate, symmetric
      tsCount1 <- floor(runif(1, min = 2, max = 20))
      tsCount2 <- tsCount1
      tsList1 <- lapply(1:tsCount1, function(x) rnorm(tsLength))
      tsList2 <- tsList1
      distMatComputed <- compDistTSList(tsList1, symbolCount = symbolCount)
    } else if (i %% 4 == 2) {# multi-variate, symmetric
      tsCount1 <- floor(runif(1, min = 2, max = 10))
      tsCount2 <- tsCount1
      tsDim <- floor(runif(1, min = 2, max = 5))
      tsList1 <- lapply(1:tsCount1, function(x) matrix(rnorm(tsLength * tsDim), ncol = tsDim))
      tsList2 <- tsList1
      distMatComputed <- compDistTSList(tsList1, symbolCount = symbolCount)
    } else if (i %% 4 == 3) {# uni-variate, asymmetric
      tsCount1 <- floor(runif(1, min = 2, max = 20))
      tsCount2 <- floor(runif(1, min = 2, max = 20))
      tsList1 <- lapply(1:tsCount1, function(x) rnorm(tsLength))
      tsList2 <- lapply(1:tsCount2, function(x) rnorm(tsLength))
      distMatComputed <- compDistTSList(list(tsList1, tsList2), symbolCount = symbolCount)
    } else if (i %% 4 == 0) {# multi-variate, asymmetric
      tsCount1 <- floor(runif(1, min = 2, max = 10))
      tsCount2 <- floor(runif(1, min = 2, max = 10))
      tsDim <- floor(runif(1, min = 2, max = 5))
      tsList1 <- lapply(1:tsCount1, function(x) matrix(rnorm(tsLength * tsDim), ncol = tsDim))
      tsList2 <- lapply(1:tsCount2, function(x) matrix(rnorm(tsLength * tsDim), ncol = tsDim))
      distMatComputed <- compDistTSList(list(tsList1, tsList2), symbolCount = symbolCount)
    }
    saxLimits <- SAXLimitsDataAdaptive(tsList = tsList2, alphabetSize = symbolCount)
    distMatExpected <- matrix(NA, nrow = tsCount1, ncol = tsCount2)
    for (k in 1:tsCount1) {
      for (l in 1:tsCount2) {
        distMatExpected[k,l] <- compDist(tsList1[[k]], tsList2[[l]], symbolLimits = saxLimits)
      }
    }
    expect_equal(distMatComputed, distMatExpected)
  }
})

test_that("pdcEntropyHeuristic same as pdc::entropyHeuristic", {
  for (i in 1:50) {
    tsCount <- floor(runif(1, min = 2, max = 20))
    # pdc::pdcDist has problems with very short ts in the entropy heuristic
    tsLength <- floor(runif(1, min = 10, max = 50))
    tsList <- lapply(1:tsCount, function(i) rnorm(tsLength))
    tsMatrix <- cbind(sapply(tsList, identity))
    # pdc::pdcDist does no check length 2 subsequences in its heuristic
    expect_equal(pdcEntropyHeuristic(tsList),
                 pdc::entropyHeuristic(X = tsMatrix, m.min = 2, m.max = 7)$m)
  }
})

test_that("pdcEntropyHeuristic works for very short time series", {
  for (i in 1:50) {
    tsCount <- floor(runif(1, min = 2, max = 20))
    expect_equal(pdcEntropyHeuristic(lapply(1:tsCount, function(i) rnorm(2))), 2)
    expect_equal(pdcEntropyHeuristic(lapply(1:tsCount, function(i) rnorm(3))), 2)
    tsLength <- floor(runif(1, min = 4, max = 10))
    m <- pdcEntropyHeuristic(lapply(1:tsCount, function(i) rnorm(tsLength)))
    expect_true(m %in% 2:min(7, tsLength - 1))
  }
})

test_that("pdcDistTwoTS same as pdc::pdcDist", {
  for (i in 1:50) {
    # pdc::pdcDist has problems with very short ts in the entropy heuristic
    tsLength <- floor(runif(1, min = 10, max = 100))
    ts1 <- rnorm(n = tsLength)
    ts2 <- rnorm(n = tsLength)
    if (i %% 2 == 1) {
      subSeqLength <- NULL
    } else {
      subSeqLength <- floor(runif(1, min = 2, max = 8))
    }
    # pdc::pdcDist does not check length 2 subsequences in its heuristic if called from pdcDist
    if (!is.null(subSeqLength) ||
        pdc::entropyHeuristic(cbind(ts1, ts2), m.min = 2, m.max = min(7, tsLength - 1))$m != 2) {
      expect_equal(pdcDistTwoTS(ts1, ts2, subSeqLength),
                   as.numeric(pdc::pdcDist(X = cbind(ts1, ts2), m = subSeqLength, t = 1)))
    }
  }
})

test_that("pdcDistTwoTSMult same as L2 norm of pdc::pdcDist to each attribute", {
  for (i in 1:50) {
    # pdc::pdcDist has problems with very short ts in the entropy heuristic
    tsLength <- floor(runif(1, min = 10, max = 50))
    tsDim <- floor(runif(1, min = 2, max = 5))
    ts1 <- matrix(rnorm(n = tsLength * tsDim), nrow = tsLength)
    ts2 <- matrix(rnorm(n = tsLength * tsDim), nrow = tsLength)
    if (i %% 2 == 1) {
      subSeqLength1 <- NULL
      # We need to execute the heuristic beforehand so it choose a global values
      # for all attributes (as in our approach), not different values for different
      # attributes
      subSeqLength2 <- pdc::entropyHeuristic(cbind(ts1, ts2), m.min = 2,
                                             m.max = min(7, tsLength - 1))$m
    } else {
      subSeqLength1 <- floor(runif(1, min = 2, max = 8))
      subSeqLength2 <- subSeqLength1
    }
    expect_equal(pdcDistTwoTSMult(ts1, ts2, subSeqLength1),
                 sqrt(sum(sapply(seq_len(tsDim), function(j) {
                   as.numeric(pdc::pdcDist(X = cbind(ts1[, j], ts2[, j]), m = subSeqLength2, t = 1))
                 })^2)))
  }
})

test_that("pdcDistList same as pdc::pdcDist", {
  for (i in 1:16) {
    tsLength <- floor(runif(1, min = 10, max = 50))
    subSeqLength <- min(floor(runif(1, min = 2, max = 8)), tsLength - 1)
    if (i %% 4 == 1) {# uni-variate, symmetric
      tsCount1 <- floor(runif(1, min = 2, max = 20))
      tsCount2 <- tsCount1
      tsList1 <- lapply(1:tsCount1, function(x) rnorm(tsLength))
      tsList2 <- tsList1
      pDistFunc <- pdcDistTwoTS
      distMatComputed <- pdcDistTSList(tsList1, subSeqLength)
    } else if (i %% 4 == 2) {# multi-variate, symmetric
      tsCount1 <- floor(runif(1, min = 2, max = 10))
      tsCount2 <- tsCount1
      tsDim <- floor(runif(1, min = 2, max = 5))
      tsList1 <- lapply(1:tsCount1, function(x) matrix(rnorm(tsLength * tsDim), ncol = tsDim))
      tsList2 <- tsList1
      pDistFunc <- pdcDistTwoTSMult
      distMatComputed <- pdcDistTSListMult(tsList1, subSeqLength)
    } else if (i %% 4 == 3) {# uni-variate, asymmetric
      tsCount1 <- floor(runif(1, min = 2, max = 20))
      tsCount2 <- floor(runif(1, min = 2, max = 20))
      tsList1 <- lapply(1:tsCount1, function(x) rnorm(tsLength))
      tsList2 <- lapply(1:tsCount2, function(x) rnorm(tsLength))
      pDistFunc <- pdcDistTwoTS
      distMatComputed <- pdcDistTSList(list(tsList1, tsList2), subSeqLength)
    } else if (i %% 4 == 0) {# multi-variate, asymmetric
      tsCount1 <- floor(runif(1, min = 2, max = 10))
      tsCount2 <- floor(runif(1, min = 2, max = 10))
      tsDim <- floor(runif(1, min = 2, max = 5))
      tsList1 <- lapply(1:tsCount1, function(x) matrix(rnorm(tsLength * tsDim), ncol = tsDim))
      tsList2 <- lapply(1:tsCount2, function(x) matrix(rnorm(tsLength * tsDim), ncol = tsDim))
      pDistFunc <- pdcDistTwoTSMult
      distMatComputed <- pdcDistTSListMult(list(tsList1, tsList2), subSeqLength)
    }
    distMatExpected <- matrix(NA, nrow = tsCount1, ncol = tsCount2)
    for (k in 1:tsCount1) {
      for (l in 1:tsCount2) {
        distMatExpected[k,l] <- do.call(pDistFunc, list(tsList1[[k]], tsList2[[l]], subSeqLength))
      }
    }
    expect_equal(distMatComputed, distMatExpected)
  }
})
