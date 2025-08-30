#' Construct Normalized Gaussian Kernel Weight Matrix
#'
#' Constructs a row-stochastic weight matrix for Nadaraya-Watson regression
#' using Gaussian kernels with specified bandwidth.
#'
#' @param x Numeric vector of predictor values.
#' @param h Numeric scalar, bandwidth parameter for the Gaussian kernel.
#'
#' @return A matrix W where W[i,j] represents the weight given to observation j
#'   when predicting at point x[i]. Each row sums to 1 (row-stochastic property).
#'
#' @details
#' The function computes a Gaussian kernel weight matrix where:
#' \deqn{K(x_i, x_j) = \frac{1}{\sqrt{2\pi}} \exp\left(-\frac{(x_i - x_j)^2}{2h^2}\right)}
#'
#' Each row is then normalized so that the weights sum to 1, making the matrix
#' row-stochastic. This ensures that the Nadaraya-Watson estimator is a proper
#' weighted average.
#'
#' @examples
#' # Simple example
#' x <- c(1, 2, 3, 4, 5)
#' h <- 0.5
#' W <- construct_W(x, h)
#' 
#' # Check row sums (should all be 1)
#' rowSums(W)
#' 
#' # Visualize weight matrix
#' image(W, main="Gaussian Kernel Weight Matrix")
#'
#' @export
construct_W <- function(x, h) {
  
  # Input validation
  if (any(is.na(x))) {
    stop("x cannot contain NA values")
  }
  if (length(h) != 1 || h <= 0) {
    stop("h must be a positive scalar")
  }
  
  n <- length(x)
  
  # Compute pairwise differences
  dist_mat <- outer(x, x, "-")
  
  # Apply Gaussian kernel
  # Using 0.3989423 ≈ 1/sqrt(2π) for computational efficiency
  K_mat <- 0.3989423 * exp(-0.5 * (dist_mat / h)^2)
  
  # Normalize each row so rows sum to 1 (row-stochastic property)
  row_sums <- rowSums(K_mat)
  
  # Avoid division by zero (though this should rarely happen with Gaussian kernels)
  if (any(row_sums == 0)) {
    warning("Some rows have zero sum. This may indicate bandwidth is too small.")
    row_sums[row_sums == 0] <- 1
  }
  
  # Return normalized matrix
  K_mat / row_sums
}