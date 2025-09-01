test_that("SNP function works correctly with large datasets", {
  # Generate test data (designed for n > 500)
  set.seed(123)
  n <- 1000
  x <- sort(runif(n, 0, 10))
  y <- sin(x) + rnorm(n, 0, 0.1)
  
  # Test basic functionality
  result <- SNP(x, y)
  
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

test_that("SNP handles different large dataset sizes", {
  sizes <- c(500, 1000, 2000)
  
  for (n in sizes) {
    set.seed(42)
    x <- sort(runif(n, 0, 5))
    y <- x^2 + rnorm(n, 0, 0.2)
    
    result <- SNP(x, y)
    
    expect_length(result$y_k_opt, n)
    expect_true(result$h_start > 0)
    expect_true(result$k_opt >= 1)
    expect_true(result$k_opt <= 10)
  }
})

test_that("SNP handles edge cases", {
  # Test input validation
  expect_error(SNP(c(1, 2), c(1)), "same length")
  expect_error(SNP(c(1, NA), c(1, 2)), "NA values")
  expect_error(SNP(c(1, 2), c(1, NA)), "NA values")
})

test_that("SNP produces reasonable smoothing on large data", {
  # Test with known smooth function on large dataset
  set.seed(456)
  n <- 800
  x <- seq(0, 2*pi, length.out = n)
  true_y <- cos(x)
  y <- true_y + rnorm(n, 0, 0.1)
  
  result <- SNP(x, y)
  
  # Smoothed result should be closer to true function than noisy observations
  mse_original <- mean((y - true_y)^2)
  mse_smoothed <- mean((result$y_k_opt - true_y)^2)
  
  expect_lt(mse_smoothed, mse_original)
})

test_that("SNP performance scales with large datasets", {
  # Test scalability - larger datasets should still run reasonably fast
  set.seed(789)
  
  times <- numeric(3)
  sizes <- c(500, 1000, 1500)
  
  for (i in seq_along(sizes)) {
    n <- sizes[i]
    x <- sort(runif(n, 0, 1))
    y <- sin(4*pi*x) + rnorm(n, 0, 0.15)
    
    start_time <- proc.time()
    result <- SNP(x, y)
    times[i] <- (proc.time() - start_time)["elapsed"]
    
    expect_length(result$y_k_opt, n)
  }
  
  # Should scale reasonably (not exponentially)
  expect_lt(times[3], times[1] * 10)  # 3x data shouldn't take 10x time
})

test_that("SNP works with heteroscedastic noise on large data", {
  # Test with varying noise levels (realistic scenario)
  set.seed(100)
  n <- 1200
  x <- sort(runif(n, 0, 8))
  
  # Create varying noise
  noise_sd <- 0.1 + 0.3 * (x > 4)  # Higher noise in second half
  y <- 2 * exp(-0.5 * x) + rnorm(n, 0, noise_sd)
  
  result <- SNP(x, y)
  
  expect_length(result$y_k_opt, n)
  expect_true(all(is.finite(result$y_k_opt)))
  expect_true(result$time_elapsed > 0)
})

test_that("SNP handles polynomial trends on large datasets", {
  # Test with polynomial data (common real-world scenario)
  set.seed(200)
  n <- 600
  x <- sort(runif(n, -2, 2))
  y <- 0.5 * x^3 - 2 * x^2 + x + 1 + rnorm(n, 0, 0.2)
  
  result <- SNP(x, y)
  
  # Should produce smooth curve
  expect_length(result$y_k_opt, n)
  expect_true(all(is.finite(result$y_k_opt)))
  
  # Smoothed curve should be less variable than original
  expect_lt(var(result$y_k_opt), var(y))
})