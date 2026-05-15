# tests/testthat/test-basic.R

test_that("construct_W works", {

  X <- matrix(runif(20), ncol = 2)

  W <- construct_W(X, h = c(0.2, 0.2))

  expect_true(is.matrix(W))

  expect_equal(nrow(W), 10)

  expect_equal(ncol(W), 10)
})


test_that("nw_direct_gcv runs", {

  x <- runif(100)

  y <- sin(2*pi*x) + rnorm(100, sd = 0.1)

  fit <- nw_direct_gcv(x, y)

  expect_true(is.list(fit))

  expect_true("y_train_opt" %in% names(fit))
})


test_that("nw_snp runs", {

  x <- runif(100)

  y <- sin(2*pi*x) + rnorm(100, sd = 0.1)

  fit <- nw_snp(x, y)

  expect_true(is.list(fit))

  expect_true("y_k_opt" %in% names(fit))
})
