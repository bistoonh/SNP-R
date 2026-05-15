# =========================================================
# dgcv.R
# =========================================================

#' Nadaraya-Watson Regression with Direct GCV Bandwidth Selection
#'
#' Performs Nadaraya-Watson kernel regression with bandwidth
#' selected by direct Generalized Cross-Validation (GCV).
#'
#' @param X Numeric matrix of predictors (n × d).
#' @param y Numeric response vector of length n.
#' @param num_h_points Number of bandwidth candidates.
#' @param mode "random" or "grid".
#' @param X_new Optional matrix of new prediction points.
#'
#' @return A list containing
#' \itemize{
#' \item y_train_opt : fitted values at training points
#' \item y_new_pred : predictions at X_new (if provided)
#' \item h_opt_gcv : optimal bandwidth vector
#' \item gcv_h : GCV values for all bandwidth candidates
#' \item h_grid : evaluated bandwidth grid
#' \item time_elapsed : computation time in seconds
#' }
#'
#' @export
nw_direct_gcv <- function(X, y, num_h_points = 30, mode = "random", X_new = NULL) {

  start_time <- proc.time()[3]

  X <- as.matrix(X)
  y <- as.numeric(y)

  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
    mode <- "grid"
  }

  if (!is.null(X_new)) {
    X_new <- as.matrix(X_new)
  }

  n <- nrow(X)
  d <- ncol(X)

  if (length(y) != n) {
    stop("Mismatch between X rows and y length")
  }

  y <- matrix(y, ncol = 1)

  # Silverman rule bandwidth range
  sd_vec <- apply(X, 2, sd)

  h_s <- 1.06 * sd_vec * n^(-1/(4 + d))

  h_min <- 0.01 * h_s
  h_max <- 2.0 * h_s

  # Generate bandwidth candidates
  if (mode == "random") {

    H_candidates <- matrix(0, nrow = num_h_points, ncol = d)

    for (j in seq_len(d)) {
      H_candidates[, j] <- runif(
        num_h_points,
        min = h_min[j],
        max = h_max[j]
      )
    }

  } else {

    grid_per_dim <- lapply(
      seq_len(d),
      function(j) seq(h_min[j], h_max[j], length.out = num_h_points)
    )

    if (d == 1) {

      H_candidates <- matrix(grid_per_dim[[1]], ncol = 1)

    } else {

      grid <- expand.grid(grid_per_dim)

      H_candidates <- as.matrix(grid)
    }
  }

  h_grid <- H_candidates
  num_combinations <- nrow(h_grid)

  gcv_h <- numeric(num_combinations)

  yhat_list <- vector("list", num_combinations)

  for (i in seq_len(num_combinations)) {

    h_current <- h_grid[i, ]

    W <- construct_W(X, h_current)

    yhat <- W %*% y

    trW <- sum(diag(W))

    rss <- sum((y - yhat)^2)

    denom <- (1 - trW / n)^2

    if (denom < 1e-12) {
      denom <- 1e-12
    }

    gcv_h[i] <- rss / denom

    yhat_list[[i]] <- yhat
  }

  inds_min <- which.min(gcv_h)

  h_opt_gcv <- h_grid[inds_min, ]

  y_train_opt <- yhat_list[[inds_min]]

  elapsed <- proc.time()[3] - start_time

  y_new <- NULL

  if (!is.null(X_new)) {

    W_new <- construct_W(X, h_opt_gcv, X_new)

    y_new <- W_new %*% y
  }

  gc()

  return(list(
    y_train_opt = y_train_opt,
    y_new_pred = y_new,
    h_opt_gcv = h_opt_gcv,
    gcv_h = gcv_h,
    h_grid = h_grid,
    time_elapsed = elapsed
  ))
}
