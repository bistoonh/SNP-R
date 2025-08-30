#' Direct Generalized Cross-Validation for Nadaraya-Watson Regression
#'
#' Implements Direct Generalized Cross-Validation (DGCV) for bandwidth selection
#' in Nadaraya-Watson regression with Gaussian kernels. This is the traditional
#' reference method that SNP aims to approximate efficiently.
#'
#' @param x Numeric vector of predictor values (sorted).
#' @param y Numeric vector of response values corresponding to x.
#' @param num_h_points Integer, number of bandwidth candidates to evaluate
#'   across the continuous bandwidth space (default: 50).
#'
#' @return A list containing:
#' \describe{
#'   \item{y_h_opt}{Final smoothed output using optimal bandwidth}
#'   \item{h_opt_gcv}{Optimal bandwidth selected by GCV}
#'   \item{gcv_h}{GCV scores for all bandwidth candidates}
#'   \item{time_elapsed}{Execution time in seconds}
#' }
#'
#' @details
#' DGCV performs an exhaustive search over a continuous bandwidth space,
#' evaluating the GCV criterion for each candidate bandwidth. While statistically
#' rigorous, this approach becomes computationally prohibitive for large datasets
#' due to the need to construct and evaluate the full smoothing matrix for each
#' bandwidth candidate.
#'
#' The bandwidth search range is determined using Silverman's rule of thumb,
#' with candidates uniformly distributed between 0.1% and 100% of the
#' Silverman bandwidth.
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
#' # Apply DGCV smoothing
#' result <- DGCV(x, y)
#' plot(x, y, pch=16, col="gray")
#' lines(x, result$y_h_opt, col="blue", lwd=2)
#'
#' # Compare with SNP
#' snp_result <- SNP(x, y)
#' lines(x, snp_result$y_k_opt, col="red", lwd=2)
#' legend("topright", c("DGCV", "SNP"), col=c("blue", "red"), lwd=2)
#'
#' @export
DGCV <- function(x, y, num_h_points = 50) {
  
  start_time <- proc.time()  # Record start time of function execution
  
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
  
  n <- length(x)
  
  # Bandwidth range based on Silverman's rule
  h_s <- 1.06 * stats::sd(x) * n^(-1/5)
  h_min <- 0.001 * h_s
  h_max <- 1 * h_s
  h_candidates <- seq(h_min, h_max, length.out = num_h_points)
  
  cat(sprintf("h_candidates: [%.4f , %.4f]\n", min(h_candidates), max(h_candidates)))
  
  # Initialize containers for results
  yk_list <- list()        # Store smoothed outputs for each h
  gcv_h <- numeric()       # Store GCV scores for each h
  
  # Loop over all bandwidth candidates
  for (i in 1:length(h_candidates)) {
    W <- construct_W(x, h_candidates[i])          # Construct Gaussian kernel weight matrix
    yhat <- W %*% y                               # Apply smoothing
    gcv_h[i] <- sum((y - yhat)^2) / ((1 - sum(diag(W)) / n)^2)  # GCV score
    yk_list[[i]] <- as.numeric(yhat)              # Save smoothed result
  }
  
  # Select h that minimizes GCV
  inds_min <- which.min(gcv_h)
  h_opt_gcv <- h_candidates[inds_min]
  
  elapsed <- proc.time() - start_time
  
  # Print summary for the user
  cat("\n--- Original GCV Smoothing Summary ---\n")
  cat("h_opt_gcv:", h_opt_gcv, "\n")
  cat("time_elapsed:", elapsed["elapsed"], "\n")
  
  gc()  # Trigger garbage collection
  
  # Return the best result and associated values
  return(list(
    y_h_opt = as.vector(yk_list[[inds_min]]),  # Final smoothed output
    h_opt_gcv = h_opt_gcv,                     # Optimal bandwidth
    gcv_h = gcv_h,                             # All GCV scores
    time_elapsed = elapsed["elapsed"]          # Elapsed time in seconds
  ))
}