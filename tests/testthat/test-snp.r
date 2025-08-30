test_that("SNP function works correctly", {
  # Generate test data
  set.seed(123)
  n <- 50
  x <- sort(runif(n, 0, 1))
  y <- sin(2*pi*x) + rnorm(n, 0, 0.1)
  
  # Test basic functionality
  result <- SNP(x, y, num_h_points = 10)  # Use fewer points for speed
  
  # Check return structure
  expect_type(result, "list")
  expect_named(result, c("y_k_opt", "h_start", "k_opt", "gcv_approx_k", "time_elapsed"))
  
  # Check dimensions
  expect_length(result$y_k_opt, n)
  expect_length(result$h_start, 1)
  expect_length(result$k_opt, 1)
  expect_type(result$gcv_approx_k, "double")
  expect_length(result$time_elapsed, 1)
  
  # Check value ranges
  expect_true(result$h_start > 0)
  expect_true(result$k_opt >= 1)
  expect_true(result$k_opt <= 10)  # k_max is 10
  expect_true(result$time_elapsed >= 0)
  expect_true(all(is.finite(result$y_k_opt)))
  expect_true(all(is.finite(result$gcv_approx_k)))
})

test_that("SNP handles edge cases", {
  # Test with minimal data
  x <- c(1, 2, 3)
  y <- c(1, 4, 2)
  
  expect_no_error(result <- SNP(x, y, num_h_points = 5))
  expect_length(result$y_k_opt, 3)
  
  # Test input validation
  expect_error(SNP(c(1, 2), c(1)), "same length")
  expect_error(SNP(c(1, NA), c(1, 2)), "NA values")
  expect_error(SNP(c(1, 2), c(1, NA)), "NA values")
  expect_error(SNP(c(1, 2), c(1, 2), num_h_points = 0), "positive")
})

test_that("SNP produces reasonable smoothing", {
  # Test with known smooth function
  set.seed(456)
  x <- seq(0, 1, length.out = 30)
  true_y <- sin(2*pi*x)
  y <- true_y + rnorm(30, 0, 0.05)
  
  result <- SNP(x, y, num_h_points = 10)
  
  # Smoothed result should be closer to true function than noisy observations
  mse_original <- mean((y - true_y)^2)
  mse_smoothed <- mean((result$y_k_opt - true_y)^2)
  
  expect_lt(mse_smoothed, mse_original)
})

test_that("SNP is deterministic given same inputs", {
  set.seed(789)
  x <- sort(runif(20, 0, 1))
  y <- x^2 + rnorm(20, 0, 0.1)
  
  # Run twice with same seed
  set.seed(100)
  result1 <- SNP(x, y, num_h_points = 5)
  
  set.seed(100)
  result2 <- SNP(x, y, num_h_points = 5)
  
  expect_equal(result1$y_k_opt, result2$y_k_opt)
  expect_equal(result1$h_start, result2$h_start)
  expect_equal(result1$k_opt, result2$k_opt)
})