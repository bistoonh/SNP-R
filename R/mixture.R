# =========================================================
# mixture.R (base R graphics with improved layout)
# =========================================================

# ---------------------------------------------------------
# Surface definition
# ---------------------------------------------------------

mixture_local_surface <- function(x1, x2) {
  ridge <- 1.2 * exp(-100 * (x2 - 0.4)^2)
  bump <- 1.0 * exp(-80 * ((x1 - 0.75)^2 + (x2 - 0.75)^2))
  hill <- 0.8 * exp(-6 * ((x1 - 0.30)^2 + (x2 - 0.25)^2))
  trough <- -0.7 * exp(-60 * ((x1 - 0.55)^2 + (x2 - 0.40)^2))
  ridge + bump + hill + trough
}


# ---------------------------------------------------------
# Data generator
# ---------------------------------------------------------

generate_mixture_data <- function(n, noise_sd = 0.5) {
  X1 <- runif(n)
  X2 <- runif(n)
  Ytrue <- mixture_local_surface(X1, X2)
  Y <- Ytrue + rnorm(n, sd = noise_sd)
  X <- cbind(X1, X2)
  list(X = X, Y = Y, Ytrue = Ytrue)
}


# ---------------------------------------------------------
# Main experiment
# ---------------------------------------------------------

#' Mixture Regression Experiment
#'
#' @param n Sample size
#' @param noise_sd Noise standard deviation
#' @param seed Random seed
#' @param grid_n Grid resolution
#'
#' @return list
#'
#' @export
mixture_experiment <- function(
    n = 5000,
    noise_sd = 0.5,
    seed = 111,
    grid_n = 70
) {

  # -------------------------------------------------------
  # Generate data
  # -------------------------------------------------------

  set.seed(seed)
  dat <- generate_mixture_data(n = n, noise_sd = noise_sd)
  X <- dat$X
  Y <- dat$Y
  Ytrue <- dat$Ytrue

  # -------------------------------------------------------
  # Fit models
  # -------------------------------------------------------

  set.seed(seed)
  res_dgcv <- nw_direct_gcv(X, Y, num_h_points = 30)

  set.seed(seed)
  res_snp <- nw_snp(X, Y, num_h_points = 30, num_slices = 50)

  yhat_dgcv <- as.numeric(res_dgcv$y_train_opt)
  yhat_snp <- as.numeric(res_snp$y_k_opt)

  # -------------------------------------------------------
  # Evaluation grid
  # -------------------------------------------------------

  g <- seq(0, 1, length.out = grid_n)
  grid <- expand.grid(X1 = g, X2 = g)
  Xg <- as.matrix(grid)

  X1g <- matrix(grid$X1, grid_n, grid_n)
  X2g <- matrix(grid$X2, grid_n, grid_n)

  Ztrue <- matrix(
    mixture_local_surface(Xg[,1], Xg[,2]),
    grid_n, grid_n
  )

  # -------------------------------------------------------
  # DGCV surface
  # -------------------------------------------------------

  h_opt <- res_dgcv$h_opt_gcv
  W_grid_dg <- construct_W(X, h_opt, Xg)
  Zdg <- matrix(W_grid_dg %*% Y, grid_n, grid_n)

  # -------------------------------------------------------
  # SNP surface
  # -------------------------------------------------------

  h_start <- res_snp$h_start
  Yk <- res_snp$y_k_minus_1_opt
  W_grid_snp <- construct_W(X, h_start, Xg)
  Zsnp <- matrix(W_grid_snp %*% Yk, grid_n, grid_n)

  # -------------------------------------------------------
  # Metrics
  # -------------------------------------------------------

  rmse_dg <- rmse(yhat_dgcv, Ytrue)
  rmse_sn <- rmse(yhat_snp, Ytrue)
  mape_dg <- mape_shift(Ytrue, yhat_dgcv)
  mape_sn <- mape_shift(Ytrue, yhat_snp)

  cat("DGCV time:", res_dgcv$time_elapsed, "\n")
  cat("DGCV RMSE:", rmse_dg, "\n")
  cat("DGCV MAPE:", mape_dg, "\n")
  cat("DGCV h_opt_gcv:\n")
  print(res_dgcv$h_opt_gcv)
  cat("\n")
  cat("SNP time:", res_snp$time_elapsed, "\n")
  cat("SNP RMSE:", rmse_sn, "\n")
  cat("SNP MAPE:", mape_sn, "\n")
  cat("SNP h_start:\n")
  print(res_snp$h_start)
  cat("SNP k_opt:", res_snp$k_opt, "\n")
  cat("SNP B:", res_snp$B, "\n")

  # -------------------------------------------------------
  # Create 2x2 layout with base R persp (improved settings)
  # -------------------------------------------------------

  # Set up color palette
  jet.colors <- colorRampPalette(c("blue", "cyan", "yellow", "red"))
  
  # Set up 2x2 layout
  par(mfrow = c(2, 2), mar = c(2, 2, 3, 1))
  
  # Common perspective settings for better view
  theta <- 30
  phi <- 30
  expand <- 0.6
  shade <- 0.4
  
  # 1. Simulated Data (scatter approximation with persp)
  # Create a binned surface for visualization
  breaks_x <- seq(0, 1, length.out = 30)
  breaks_y <- seq(0, 1, length.out = 30)
  
  x_bins <- cut(X[,1], breaks = breaks_x, labels = FALSE, include.lowest = TRUE)
  y_bins <- cut(X[,2], breaks = breaks_y, labels = FALSE, include.lowest = TRUE)
  
  Z_binned <- matrix(NA, length(breaks_x)-1, length(breaks_y)-1)
  for(i in 1:(length(breaks_x)-1)) {
    for(j in 1:(length(breaks_y)-1)) {
      idx <- which(x_bins == i & y_bins == j)
      if(length(idx) > 0) {
        Z_binned[i, j] <- mean(Y[idx])
      }
    }
  }
  
  # Interpolate missing values
  for(i in 1:nrow(Z_binned)) {
    for(j in 1:ncol(Z_binned)) {
      if(is.na(Z_binned[i,j])) {
        neighbors <- c()
        if(i > 1 && !is.na(Z_binned[i-1,j])) neighbors <- c(neighbors, Z_binned[i-1,j])
        if(i < nrow(Z_binned) && !is.na(Z_binned[i+1,j])) neighbors <- c(neighbors, Z_binned[i+1,j])
        if(j > 1 && !is.na(Z_binned[i,j-1])) neighbors <- c(neighbors, Z_binned[i,j-1])
        if(j < ncol(Z_binned) && !is.na(Z_binned[i,j+1])) neighbors <- c(neighbors, Z_binned[i,j+1])
        if(length(neighbors) > 0) Z_binned[i,j] <- mean(neighbors)
      }
    }
  }
  
  # Fill remaining NAs with global mean
  Z_binned[is.na(Z_binned)] <- mean(Y, na.rm = TRUE)
  
  nrz <- nrow(Z_binned)
  ncz <- ncol(Z_binned)
  nbcol <- 100
  color <- jet.colors(nbcol)
  zfacet <- Z_binned[-1, -1] + Z_binned[-1, -ncz] + Z_binned[-nrz, -1] + Z_binned[-nrz, -ncz]
  facetcol <- cut(zfacet, nbcol)
  
  persp(breaks_x[-length(breaks_x)], breaks_y[-length(breaks_y)], Z_binned,
        theta = theta, phi = phi, expand = expand, col = color[facetcol],
        shade = shade, border = NA,
        xlab = "X1", ylab = "X2", zlab = "Y",
        main = "Simulated Data\n(Mixture of Local Features)")
  
  # 2. True Surface
  nrz <- nrow(Ztrue)
  ncz <- ncol(Ztrue)
  zfacet <- Ztrue[-1, -1] + Ztrue[-1, -ncz] + Ztrue[-nrz, -1] + Ztrue[-nrz, -ncz]
  facetcol <- cut(zfacet, nbcol)
  
  persp(g, g, Ztrue,
        theta = theta, phi = phi, expand = expand, col = color[facetcol],
        shade = shade, border = NA,
        xlab = "X1", ylab = "X2", zlab = "Y",
        main = "True Surface\n(Mixture of Local Features)")
  
  # 3. DGCV Fitted Surface
  nrz <- nrow(Zdg)
  ncz <- ncol(Zdg)
  zfacet <- Zdg[-1, -1] + Zdg[-1, -ncz] + Zdg[-nrz, -1] + Zdg[-nrz, -ncz]
  facetcol <- cut(zfacet, nbcol)
  
  persp(g, g, Zdg,
        theta = theta, phi = phi, expand = expand, col = color[facetcol],
        shade = shade, border = NA,
        xlab = "X1", ylab = "X2", zlab = "Y",
        main = sprintf("DGCV Fitted Surface\nRMSE=%.4f", rmse_dg))
  
  # 4. SNP Fitted Surface
  nrz <- nrow(Zsnp)
  ncz <- ncol(Zsnp)
  zfacet <- Zsnp[-1, -1] + Zsnp[-1, -ncz] + Zsnp[-nrz, -1] + Zsnp[-nrz, -ncz]
  facetcol <- cut(zfacet, nbcol)
  
  persp(g, g, Zsnp,
        theta = theta, phi = phi, expand = expand, col = color[facetcol],
        shade = shade, border = NA,
        xlab = "X1", ylab = "X2", zlab = "Y",
        main = sprintf("SNP Fitted Surface\nRMSE=%.4f", rmse_sn))
  
  par(mfrow = c(1, 1))

  # -------------------------------------------------------
  # GCV profiles (1x2 layout)
  # -------------------------------------------------------

  gcv_k <- res_snp$gcv_approx_k
  k_vals <- seq_along(gcv_k)

  h_raw <- as.matrix(res_dgcv$h_grid)
  gcv_raw <- res_dgcv$gcv_h
  h_scalar <- sqrt(h_raw[,1] * h_raw[,2])
  idx <- order(h_scalar)
  h_sorted <- h_scalar[idx]
  gcv_sorted <- gcv_raw[idx]

  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  
  plot(k_vals, gcv_k, type = "b", col = "blue", pch = 19, lwd = 2,
       main = "GCV(k) for SNP", xlab = "k", ylab = "GCV(k)",
       cex.main = 1.2, cex.lab = 1.1)
  grid()
  
  plot(h_sorted, gcv_sorted, type = "b", col = "red", pch = 19, lwd = 2, log = "x",
       main = "GCV(h) for DGCV", xlab = "Bandwidth (log scale)", ylab = "GCV(h)",
       cex.main = 1.2, cex.lab = 1.1)
  grid()
  
  par(mfrow = c(1, 1))

  # -------------------------------------------------------
  # Return results
  # -------------------------------------------------------

  invisible(list(
    dgcv = res_dgcv,
    snp = res_snp,
    rmse_dgcv = rmse_dg,
    rmse_snp = rmse_sn,
    mape_dgcv = mape_dg,
    mape_snp = mape_sn
  ))
}
