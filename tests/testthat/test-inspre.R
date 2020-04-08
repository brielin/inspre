D <- 3
N <- 10
p_net <- 1
dataset <- generate_dataset(D, N, p_net)
res <- cor_w_se(dataset$X)

test_that("make_weights_works", {
  weights <- make_weights(res$SE_S)
  expect_equal(dim(weights), dim(res$SE_S))
})

test_that("make_weights_max_works", {
  weights <- make_weights(res$SE_S, max_min_ratio = 1.5)
  expect_equal(max(weights)/min(weights), 1.5)
})

test_that("make_weights_na_works", {
  SE <- res$SE_S
  SE[1, 2]<- NA

  weights <- make_weights(SE)
  expect_equal(sum(is.na(weights)), 0)

  weights <- make_weights(SE, max_min_ratio = 1.5)
  expect_equal(sum(is.na(weights)), 0)
  expect_equal(max(weights)/min(weights[weights > 0]), 1.5)
})

test_that("inspre_worker_works", {
  X <- res$S_hat
  SE <- res$SE_S
  X[1, 2] <- NA
  SE[1, 2] <- NA
  W <- make_weights(SE)
  res <- inspre_worker(X, W, its = 10)
  expect_equal(sum(is.na(res$V)), 0)
  expect_equal(dim(res$V), dim(X))
  expect_equal(dim(res$U), dim(X))
})

test_that("inspre_works", {
  X <- res$S_hat
  SE <- res$SE_S
  X[1, 2] <- NA
  SE[1, 2] <- NA
  W <- make_weights(SE)
  res <- inspre(X, W, its = 10, nlambda = 10, verbose = 0)
  expect_equal(sum(is.na(res$V)), 0)
  expect_equal(dim(res$V), c(dim(X), 10))
  expect_equal(length(res$lambda), 10)
  expect_equal(length(res$rho), 10)
  expect_equal(length(res$L), 10)
})

test_that("inspre_cv_works", {
  X <- res$S_hat
  SE <- res$SE_S
  X[1, 2] <- NA
  SE[1, 2] <- NA
  W <- make_weights(SE)
  res <- inspre(X, W, its = 10, nlambda = 10, verbose = 0, cv_folds = 3)
  expect_equal(length(res$D_hat), 10)
  expect_equal(length(res$D_hat_se), 10)
})
