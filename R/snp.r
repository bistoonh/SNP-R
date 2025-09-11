#' Stepwise Noise Peeling for Nadaraya-Watson Regression
#'
#' Implements the Stepwise Noise Peeling (SNP) algorithm that bypasses 
#' bandwidth selection in Nadaraya-Watson regression by using iterative 
#' smoothing. SNP provides a scalable alternative to Direct Generalized 
#' Cross-Validation (DGCV) by avoiding continuous bandwidth optimization.
#'
#' @param x Numeric vector of predictor values (sorted).
#' @param y Numeric vector of response values corresponding to x.
#'
#' @return A list containing:
#' \describe{
#'   \item{y_k_opt}{Final smoothed output vector}
#'   \item{h_start}{Final chosen initial bandwidth}
#'   \item{k_opt}{Optimal number of iterations}
#'   \item{gcv_approx_k}{GCV values for each iteration}
#'   \item{time_elapsed}{Execution time in seconds}
#' }
#'
#' @details
#' The SNP algorithm operates in two phases:
#' \describe{
#'   \item{Phase I}{Constructs a conservative initial bandwidth using random 
#'     slices of data and lightweight GCV within each slice}
#'   \item{Phase II}{Fixes the smoothing operator and repeatedly applies it, 
#'     selecting optimal iterations via discrete GCV}
#' }
#'
#' The algorithm converts costly continuous bandwidth search into lightweight 
#' discrete selection, reducing runtime by orders of magnitude while yielding
#' estimates statistically equivalent to DGCV.
#'
#' @references
#' Your paper reference here (to be added after publication)
#'
#' @examples
#' # Generate sample data
#' set.seed(123)
#' n <- 100
#' x <- sort(runif(n, 0, 1))
#' y <- sin(2*pi*x) + rnorm(n, 0, 0.1)
#'
#' # Apply SNP smoothing
#' result <- SNP(x, y)
#' plot(x, y, pch=16, col="gray")
#' lines(x, result$y_k_opt, col="red", lwd=2)
#'
#' @export
SNP <- function(x, y) {
  start_time <- proc.time()   # Record start time of function execution
  n <- length(x)
  
  num_h_points <- 40
  num_slices <- 60
  k_max <- 10
  
  # Input validation
  if (length(x) != length(y)) {
    stop("x and y must have the same length")
  }
  if (any(is.na(x)) || any(is.na(y))) {
    stop("x and y cannot contain NA values")
  }
  if (num_h_points <= 0) {
    stop("num_h_points must be positive")
  }

  
  # Initial bandwidth range based on Silverman's rule of thumb
  h_s <- 1.06 * stats::sd(x) * n^(-1/5)
  # Lower bound for bandwidth: ensures W is not too sparse
  h_min <- 0.001 * h_s
  # Upper bound for bandwidth: standard Silverman bandwidth
  h_max <- 1 * h_s
  
  cat("-------------Start SNP-------------\n")
  cat(sprintf("h_candidates: [%.4f , %.4f]\n", h_min, h_max))
  
  # Determine slice size based on sample size
  min_slice <- 50
  
  if (n < min_slice) {
    slice_size <- n
  } else {
    #slice_size <- floor(max(min_slice, sqrt(n * log(n))))
    #slice_size <- floor(sqrt(n * log(n)))
  }
  slice_size <- floor(sqrt(n * log(n)))
  # Randomly select starting indices for each slice
  start_indices <- sample(1:(n - slice_size + 1), num_slices, replace = TRUE)
  slice_indices <- lapply(start_indices, function(start_idx) start_idx:(start_idx + slice_size - 1))
  
  # Internal function to compute optimal h for a given slice using GCV
  compute_h_opt <- function(idx) {
    x_slice <- x[idx]
    y_slice <- y[idx]
    h_candidates <- stats::runif(num_h_points, h_min, h_max)
    # Compute GCV for each h candidate
    gcv_scores <- sapply(h_candidates, function(h) {
      W_slice <- construct_W(x_slice, h)
      y_hat <- W_slice %*% y_slice
      mean((y_slice - y_hat)^2) / ((1 - mean(diag(W_slice)))^2)
    })
    
    # Return the h with minimum GCV
    h_candidates[which.min(gcv_scores)]
  }
  
  # Apply to all slices
  h_opts <- sapply(slice_indices, compute_h_opt)
  
  elapsed_h <- proc.time() - start_time
  start_time <- proc.time()
  
  # Use median of h_opt estimates as starting point
  h_start <- 0.5 * stats::median(h_opts)
  
  
  cat("h_start:", h_start, "\n")
  cat("summary h_opts:", summary(h_opts), "\n")
  
  # Trace approximation function
  trace_Wk <- function(trWh, k, cap_one = TRUE) {
    val <- 1 + (trWh - 1) / sqrt(k)
    if (cap_one) pmax(1, val) else val
  }
  
  # Adaptive h_start adjustment
  i0 <- 1
  while (i0 <= 10) {
    if (i0 == 10) cat("Last chance for change h_start\n")
    
    # Initial weight matrix with h_start
    W <- construct_W(x, h_start)
    trWh <- sum(diag(W))
    
    # Apply initial smoothing
    yk <- W %*% y
    yk_list <- list()
    gcv_approx_k <- numeric()
    traces <- numeric(k_max)
    
    for (k in 1:k_max) {
      traces[k] <- trace_Wk(trWh, k)
      gcv_approx_k[k] <- sum((y - yk)^2) / ((1 - traces[k] / n)^2)
      yk_list[[k]] <- as.numeric(yk)
      yk <- W %*% yk  # Apply one more smoothing iteration
    }
    
    # Choose k (number of iterations) with lowest approximate GCV
    k_opt <- which.min(gcv_approx_k)
    
    if (length(k_opt) == 1) {
      if (k_opt == 1) {
        h_start <- 0.5 * h_start
        i0 <- i0 + 1
        cat("new smaller h_start:", h_start, "\n")
      } else if (k_opt == k_max) {
        h_start <- 0.5 * sqrt(k_max) * h_start
        i0 <- i0 + 1
        cat("new bigger h_start:", h_start, "\n")
      } else {
        i0 <- 11
      }
    } else {
      cat("new bigger h_start:", h_start, " because trace = 0\n")
      h_start <- 1.5 * h_start
      i0 <- max(1, i0 - 2)
    }
  }
  
  elapsed_k <- proc.time() - start_time
  elapsed <- elapsed_h + elapsed_k
  
  # Print summary
  cat("\n--- Adaptive Smoothing Summary ---\n")
  cat("time_elapsed_h:", elapsed_h["elapsed"], ', time_elapsed_k:', elapsed_k["elapsed"], "\n")
  cat("h_start (final):", h_start, "\n")
  cat("k_opt (final):", k_opt, "\n")
  cat("k_max:", k_max, "\n")
  cat("time_elapsed:", elapsed["elapsed"], "\n")
  cat("\n")
  cat("-------------End SNP-------------\n")
  
  gc()  # Trigger garbage collection
  
  # Return results
  return(list(
    y_k_opt = as.vector(yk_list[[k_opt]]),  # Final smoothed output
    h_start = h_start,                      # Final chosen h
    k_opt = k_opt,                          # Optimal iteration count
    gcv_approx_k = gcv_approx_k,            # GCV values per iteration
    time_elapsed = elapsed["elapsed"]       # Elapsed time in seconds
  ))
}
