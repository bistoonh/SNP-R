# =========================================================
# snp.R
# =========================================================

# ---------------------------------------------------------
# Generate local slices using nearest neighbors
# ---------------------------------------------------------
generate_slices <- function(X, s, T) {

  X <- as.matrix(X)
  n <- nrow(X)

  slices <- vector("list", T)

  for (t in seq_len(T)) {

    center_idx <- sample.int(n, 1)

    center <- X[center_idx, ]

    distances <- sqrt(rowSums((X - matrix(center, n, ncol(X), byrow = TRUE))^2))

    slice_indices <- order(distances)[1:s]

    slices[[t]] <- slice_indices
  }

  return(slices)
}


# ---------------------------------------------------------
# SNP Nadaraya-Watson
# ---------------------------------------------------------

#' Stepwise Noise Peeling for Nadaraya-Watson Regression
#'
#' @param X Predictor matrix (n × d)
#' @param y Response vector
#' @param num_h_points Number of bandwidth candidates
#' @param num_slices Number of slices
#' @param X_new Optional prediction points
#'
#' @return list containing fitted values, bandwidth and diagnostics
#'
#' @export
nw_snp <- function(X, y, num_h_points = 30, num_slices = 50, X_new = NULL) {

  X <- as.matrix(X)
  y <- as.numeric(y)

  if (is.vector(X)) {
    X <- matrix(X, ncol = 1)
  }

  y <- matrix(y, ncol = 1)

  if (!is.null(X_new)) {
    X_new <- as.matrix(X_new)
  }

  start_time <- proc.time()[3]

  k_max <- 10
  rho <- 0.5

  n <- nrow(X)
  d <- ncol(X)

  # Silverman-like bandwidth
  sd_vec <- apply(X, 2, sd)

  h_s <- 1.06 * sd_vec * n^(-1/(4 + d))

  h_min <- 0.01 * h_s
  h_max <- 2.0 * h_s

  # Slice size
  min_slice <- 50

  if (n < min_slice) {
    slice_size <- n
  } else {
    slice_size <- floor(max(min_slice, sqrt(n * log(n))))
  }

  slice_indices <- generate_slices(X, slice_size, num_slices)

  # -----------------------------------------------------
  # Phase I : bandwidth initialization
  # -----------------------------------------------------

  compute_h_opt <- function(idx) {

    X_slice <- X[idx, , drop = FALSE]
    y_slice <- y[idx, , drop = FALSE]

    h_candidates <- matrix(0, num_h_points, d)

    for (j in seq_len(d)) {
      h_candidates[, j] <- runif(
        num_h_points,
        h_min[j],
        h_max[j]
      )
    }

    gcv_scores <- numeric(num_h_points)

    for (i in seq_len(num_h_points)) {

      h_vec <- h_candidates[i, ]

      W_slice <- construct_W(X_slice, h_vec)

      y_hat <- W_slice %*% y_slice

      denom <- (1 - mean(diag(W_slice)))^2

      if (denom < 1e-12) {
        denom <- 1e-12
      }

      gcv_scores[i] <- mean((y_slice - y_hat)^2) / denom
    }

    return(h_candidates[which.min(gcv_scores), ])
  }

  h_opts <- do.call(
    rbind,
    lapply(slice_indices, compute_h_opt)
  )

  elapsed_h <- proc.time()[3] - start_time

  h_start <- rho * apply(h_opts, 2, median)

  # -----------------------------------------------------
  # Trace approximation
  # -----------------------------------------------------

  trace_Wk <- function(trWh, k, d) {

    pmax(
      1,
      1 + (trWh - 1) / (k^(d / 2))
    )
  }

  # -----------------------------------------------------
  # Phase II : iterative smoothing
  # -----------------------------------------------------

  start_time_k <- proc.time()[3]

  max_restarts <- 10
  restart_count <- 0
  B <- 0

  while (restart_count < max_restarts) {

    B <- B + 1

    if (restart_count == max_restarts - 1) {
      message("Last chance for h_start adjustment")
    }

    W <- construct_W(X, h_start)

    trWh <- sum(diag(W))

    yk <- W %*% y

    yk_list <- vector("list", k_max)

    gcv_approx_k <- numeric(k_max)

    traces <- numeric(k_max)

    for (k in seq_len(k_max)) {

      traces[k] <- trace_Wk(trWh, k, d)

      denom <- (1 - traces[k] / n)^2

      if (denom < 1e-12) {
        denom <- 1e-12
      }

      gcv_approx_k[k] <- sum((y - yk)^2) / denom

      yk_list[[k]] <- yk

      if (k < k_max) {
        yk <- W %*% yk
      }
    }

    k_opt <- which.min(gcv_approx_k)

    if (k_opt == 1) {

      h_start <- rho * h_start

      restart_count <- restart_count + 1

      message("Restart: smaller h_start")
      print(h_start)

    } else if (k_opt == k_max) {

      h_start <- rho * h_start * sqrt(k_max)

      restart_count <- restart_count + 1

      message("Restart: larger h_start")
      print(h_start)

    } else {
      break
    }
  }

  elapsed_k <- proc.time()[3] - start_time_k

  elapsed <- elapsed_h + elapsed_k

  # -----------------------------------------------------
  # Prediction
  # -----------------------------------------------------

  y_new_pred <- NULL

  if (!is.null(X_new)) {

    W_new <- construct_W(X, h_start, X_new)

    if (k_opt == 1) {
      y_prev <- y
    } else {
      y_prev <- yk_list[[k_opt - 1]]
    }

    y_new_pred <- W_new %*% y_prev
  }

  gc()

  return(list(
    y_k_opt = yk_list[[k_opt]],
    y_k_minus_1_opt = if (k_opt == 1) y else yk_list[[k_opt - 1]],
    y_new_pred = y_new_pred,
    h_start = h_start,
    k_opt = k_opt,
    gcv_approx_k = gcv_approx_k,
    traces = traces,
    time_elapsed = elapsed,
    B = B
  ))
}
