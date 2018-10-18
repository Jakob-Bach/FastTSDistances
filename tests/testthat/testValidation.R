context("Cluster Validity Indices")

test_that("entropy_fast equivalent to routine in mcclust::vi.dist()", {
  # code copied from mcclust::vi.dist() [is a nested function there]
  ent <- function(cl) {
    n <- length(cl)
    p <- table(cl)/n
    -sum(p * log(p, base = 2))
  }
  for (i in 1:100) {
    n <- floor(runif(1, min = 20, max = 50))
    k <- floor(runif(1, min = 2, max = 10))
    # Fake cluster assignments; each number at least once
    assignments <- c(sample(1:k, size = k), sample(1:k, size = n - k, replace = TRUE))
    expect_equal(clusterEntropy_fast(assignments), ent(assignments))
  }
})

test_that("generalizedDB_fast same as clv::clv.Davies.Bouldin()", {
  for (i in 1:100) {
    n <- floor(runif(1, min = 20, max = 50))
    k <- floor(runif(1, min = 2, max = 10))
    data <- rnorm(n)
    distMat <- as.matrix(proxy::dist(data))
    # Fake cluster assignments; each number at least once
    assignments <- c(sample(1:k, size = k), sample(1:k, size = n - k, replace = TRUE))
    interIntraDistances <- clv::cls.scatt.diss.mx(diss.mx = distMat, clust = assignments)
    expect_equal(generalizedDB_fast(interIntraDistances$intercls.average,
                    interIntraDistances$intracls.average),
                 as.numeric(clv::clv.Davies.Bouldin(interIntraDistances, intracls = "average",
                    intercls = "average")))
  }
})

test_that("generalizedDunn_fast same as clv::clv.Dunn()", {
  for (i in 1:100) {
    n <- floor(runif(1, min = 20, max = 50))
    k <- floor(runif(1, min = 2, max = 10))
    data <- rnorm(n)
    distMat <- as.matrix(proxy::dist(data))
    # Fake cluster assignments; each number at least once
    assignments <- c(sample(1:k, size = k), sample(1:k, size = n - k, replace = TRUE))
    interIntraDistances <- clv::cls.scatt.diss.mx(diss.mx = distMat, clust = assignments)
    expect_equal(generalizedDunn_fast(interIntraDistances$intercls.average,
                                    interIntraDistances$intracls.average),
                 as.numeric(clv::clv.Dunn(interIntraDistances, intracls = "average",
                                                    intercls = "average")))
  }
})

test_that("conditionalEntropy_fast (normalized) in [0,1]", {
  expect_equal(conditionalEntropy_fast(c(1,1,2,2), c(1,2,1,2), normalizeAndInvert = TRUE), 0)
  expect_equal(conditionalEntropy_fast(assignments = c(1,1,2,2,3,3,4,4),
                                       groundTruth = c(1,1,1,1,2,2,2,2), normalizeAndInvert = TRUE), 1)
  for (i in 1:100) {
    n <- floor(runif(1, min = 20, max = 50))
    k <- floor(runif(1, min = 2, max = 10))
    # Fake cluster assignments; each number at least once
    assignments1 <- as.integer(c(sample(1:k, size = k), sample(1:k, size = n - k, replace = TRUE)))
    groundTruth <- as.integer(c(sample(1:k, size = k), sample(1:k, size = n - k, replace = TRUE)))
    condEntropy <- conditionalEntropy_fast(assignments1, groundTruth, normalizeAndInvert = TRUE)
    expect_true(0 <= condEntropy && condEntropy <= 1)
    expect_equal(vanDongen_fast(assignments1, assignments1, normalizeAndInvert = TRUE), 1)
  }
})

test_that("randIndex_fast equivalent to clv::clv.Rand()", {
  for (i in 1:50) {
    n <- floor(runif(1, min = 20, max = 50))
    k <- floor(runif(1, min = 2, max = 10))
    # Fake cluster assignments; each number at least once
    assignments1 <- as.integer(c(sample(1:k, size = k), sample(1:k, size = n - k, replace = TRUE)))
    assignments2 <- as.integer(c(sample(1:k, size = k), sample(1:k, size = n - k, replace = TRUE)))
    expect_equal(randIndex_fast(pairCVIParameters_fast(assignments1, assignments2)),
                 clv::clv.Rand(clv::std.ext(assignments1, assignments2)))
  }
})

test_that("randIndex_fast (normalized) equivalent to igraph::compare", {
  for (i in 1:50) {
    n <- floor(runif(1, min = 20, max = 50))
    k <- floor(runif(1, min = 2, max = 10))
    # Fake cluster assignments; each number at least once
    assignments1 <- as.integer(c(sample(1:k, size = k), sample(1:k, size = n - k, replace = TRUE)))
    assignments2 <- as.integer(c(sample(1:k, size = k), sample(1:k, size = n - k, replace = TRUE)))
    expect_equal(randIndex_fast(pairCVIParameters_fast(assignments1, assignments2), normalize = TRUE),
                 igraph::compare(assignments1, assignments2, method = "adjusted.rand"))
  }
})

test_that("fowlkesMallows_fast equivalent to clv::clv.Fowlkes.Mallows()", {
  for (i in 1:100) {
    n <- floor(runif(1, min = 20, max = 50))
    k <- floor(runif(1, min = 2, max = 10))
    # Fake cluster assignments; each number at least once
    assignments1 <- as.integer(c(sample(1:k, size = k), sample(1:k, size = n - k, replace = TRUE)))
    assignments2 <- as.integer(c(sample(1:k, size = k), sample(1:k, size = n - k, replace = TRUE)))
    expect_equal(fowlkesMallows_fast(pairCVIParameters_fast(assignments1, assignments2)),
                 clv::clv.Folkes.Mallows(clv::std.ext(assignments1, assignments2)))
  }
})

test_that("phi_fast equivalent to clusterCrit::extCriteria", {
  for (i in 1:100) {
    n <- floor(runif(1, min = 20, max = 50))
    k <- floor(runif(1, min = 2, max = 10))
    # Fake cluster assignments; each number at least once
    assignments1 <- as.integer(c(sample(1:k, size = k), sample(1:k, size = n - k, replace = TRUE)))
    assignments2 <- as.integer(c(sample(1:k, size = k), sample(1:k, size = n - k, replace = TRUE)))
    # a bit more numerical differences than in other tests, therefore manual tolerance
    # naming of coefficient varies, in clusterCrit its called Hubert
    expect_true(abs(phi_fast(pairCVIParameters_fast(assignments1, assignments2)) -
                 as.numeric(clusterCrit::extCriteria(assignments1, assignments2, crit = "Hubert"))) <= 1e-6)
  }
})

test_that("purity correct for test cases", {
  expect_equal(purity_fast(assignments = c(1,1,2,2,2), groundTruth = c(2,2,1,1,1)), 1)
  expect_equal(purity_fast(assignments = 1:5, groundTruth = c(2,2,1,1,1)), 1)
  expect_false(purity_fast(groundTruth = 1:5, assignments = c(2,2,1,1,1)) == 1)
  expect_equal(purity_fast(assignments = c(1,1,2,2,2), groundTruth = c(2,2,1,3,1)), 0.8)
  expect_equal(purity_fast(assignments = c(1,1,1,1,2,2,2,2), groundTruth = c(1,2,1,2,1,2,1,2)), 0.5)
})

test_that("vanDongen_fast equivalent to igraph::compare", {
  for (i in 1:50) {
    n <- floor(runif(1, min = 20, max = 50))
    k <- floor(runif(1, min = 2, max = 10))
    # Fake cluster assignments; each number at least once
    assignments1 <- as.integer(c(sample(1:k, size = k), sample(1:k, size = n - k, replace = TRUE)))
    assignments2 <- as.integer(c(sample(1:k, size = k), sample(1:k, size = n - k, replace = TRUE)))
    expect_equal(vanDongen_fast(assignments1, assignments2),
                 igraph::compare(assignments1, assignments2, method = "split.join"))
  }
})

test_that("vanDongen_fast (normalized) in [0,1]", {
  expect_equal(vanDongen_fast(c(1,1,2,2), c(1,2,1,2), normalizeAndInvert = TRUE), 0)
  for (i in 1:50) {
    n <- floor(runif(1, min = 20, max = 50))
    k <- floor(runif(1, min = 2, max = 10))
    # Fake cluster assignments; each number at least once
    assignments1 <- as.integer(c(sample(1:k, size = k), sample(1:k, size = n - k, replace = TRUE)))
    assignments2 <- as.integer(c(sample(1:k, size = k), sample(1:k, size = n - k, replace = TRUE)))
    vanDongen <- vanDongen_fast(assignments1, assignments2, normalizeAndInvert = TRUE)
    expect_true(0 <= vanDongen && vanDongen <= 1)
    expect_equal(vanDongen_fast(assignments1, assignments1, normalizeAndInvert = TRUE), 1)
  }
})

test_that("VI_fast equivalent to mcclust::vi.dist", {
  for (i in 1:50) {
    n <- floor(runif(1, min = 20, max = 50))
    k <- floor(runif(1, min = 2, max = 10))
    # Fake cluster assignments; each number at least once
    assignments1 <- as.integer(c(sample(1:k, size = k), sample(1:k, size = n - k, replace = TRUE)))
    assignments2 <- as.integer(c(sample(1:k, size = k), sample(1:k, size = n - k, replace = TRUE)))
    expect_equal(VI_fast(assignments1, assignments2),
                 mcclust::vi.dist(assignments1, assignments2, base = 2))
  }
})

test_that("VI_fast (normalized) equivalent to igraph::compare", {
  for (i in 1:50) {
    n <- floor(runif(1, min = 20, max = 50))
    k <- floor(runif(1, min = 2, max = 10))
    # Fake cluster assignments; each number at least once
    assignments1 <- as.integer(c(sample(1:k, size = k), sample(1:k, size = n - k, replace = TRUE)))
    assignments2 <- as.integer(c(sample(1:k, size = k), sample(1:k, size = n - k, replace = TRUE)))
    expect_equal(VI_fast(assignments1, assignments2, normalizeAndInvert = TRUE),
                 igraph::compare(assignments1, assignments2, method = "nmi"))
  }
})
