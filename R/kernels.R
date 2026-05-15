# =========================================================
# kernels.R
# =========================================================

#' Construct Nadaraya-Watson Smoother Matrix
#'
#' Computes the Nadaraya-Watson smoother matrix using a
#' Gaussian kernel.
#'
#' @param X_train Numeric matrix of training predictors
#'   with shape (n_train, d).
#' @param h Numeric scalar or vector of bandwidths.
#' @param X_new Optional numeric matrix of new predictor
#'   points with shape (n_new, d). If NULL, X_train is used.
#'
#' @return A numeric matrix of normalized kernel weights.
#'   Each row sums to 1.
#'
#' @details
#' The Gaussian kernel is
#'
#' \deqn{
#' K(u) = \exp(-0.5 ||u||^2)
#' }
#'
#' and the smoother matrix elements are
#'
#' \deqn{
#' W_{ij} =
#' \frac{
#' K((X_{new,i} - X_{train,j}) / h)
#' }{
#' \sum_k K((X_{new,i} - X_{train,k}) / h)
#' }
#' }
#'
#' @examples
#' X <- matrix(c(1, 2, 3), ncol = 1)
#'
#' W <- construct_W(X, h = 1)
#'
#' dim(W)
#' rowSums(W)
#'
#' @export
# =========================================================
# Fast Gaussian NW smoothing matrix
# =========================================================

construct_W <- function(X, h, X_new = NULL) {

  X <- as.matrix(X)

  if (is.null(X_new)) {
    X_new <- X
  } else {
    X_new <- as.matrix(X_new)
  }

  h <- as.numeric(h)

  # -------------------------------------------------------
  # Scale data by bandwidth
  # -------------------------------------------------------

  Xs <- sweep(X, 2, h, "/", check.margin = FALSE)
  Zs <- sweep(X_new, 2, h, "/", check.margin = FALSE)

  # -------------------------------------------------------
  # Squared norms
  # -------------------------------------------------------

  X_norm <- rowSums(Xs * Xs)
  Z_norm <- rowSums(Zs * Zs)

  # -------------------------------------------------------
  # Pairwise squared distances
  #
  # ||x-z||² = ||x||² + ||z||² - 2 xᵀz
  # -------------------------------------------------------

  G <- Zs %*% t(Xs)

  D2 <- matrix(
    Z_norm,
    nrow = length(Z_norm),
    ncol = length(X_norm)
  )

  D2 <- D2 + rep(X_norm, each = length(Z_norm))

  D2 <- D2 - 2 * G

  # numerical stability
  D2[D2 < 0] <- 0

  # -------------------------------------------------------
  # Gaussian kernel
  # -------------------------------------------------------

  K <- exp(-0.5 * D2)

  # -------------------------------------------------------
  # Row normalization
  # -------------------------------------------------------

  rs <- rowSums(K)

  rs[rs == 0] <- 1

  W <- K / rs

  return(W)
}

