# =========================================================
# realdata_2d.R
# =========================================================

# ---------------------------------------------------------
# Load California housing data
# ---------------------------------------------------------

load_realdata_2d <- function() {

  data("housing", package = "snpregR")

  df <- na.omit(housing)

  x1 <- df$median_income
  x2 <- df$housing_median_age
  y  <- df$median_house_value

  X <- cbind(x1, x2)

  list(
    X = X,
    x1 = x1,
    x2 = x2,
    y = y
  )
}


# ---------------------------------------------------------
# Main experiment
# ---------------------------------------------------------

#' California Housing Example (2D)
#'
#' Fits SNP-NW and Direct-GCV NW using:
#' median_income, housing_median_age -> median_house_value
#'
#' @param seed Random seed
#'
#' @return list
#'
#' @export
realdata_2d <- function(seed = 111) {

  dat <- load_realdata_2d()

  X  <- dat$X
  x1 <- dat$x1
  x2 <- dat$x2
  y  <- dat$y

  # -------------------------------------------------------
  # DGCV
  # -------------------------------------------------------

  set.seed(seed)

  res_dgcv <- nw_direct_gcv(
    X,
    y,
    num_h_points = 30,
    mode = "random"
  )

  # -------------------------------------------------------
  # SNP
  # -------------------------------------------------------

  set.seed(seed)

  res_snp <- nw_snp(
    X,
    y,
    num_h_points = 30,
    num_slices = 50
  )

  yhat_dgcv <- as.numeric(res_dgcv$y_train_opt)
  yhat_snp  <- as.numeric(res_snp$y_k_opt)

  rmse_dg <- rmse(y, yhat_dgcv)
  rmse_sn <- rmse(y, yhat_snp)

  cat("DGCV time:", res_dgcv$time_elapsed, "\n")
  cat("DGCV RMSE:", rmse_dg, "\n")
  cat("DGCV h_opt_gcv:\n")
  print(res_dgcv$h_opt_gcv)

  cat("\n")

  cat("SNP time:", res_snp$time_elapsed, "\n")
  cat("SNP RMSE:", rmse_sn, "\n")
  cat("SNP h_start:\n")
  print(res_snp$h_start)
  cat("SNP k_opt:", res_snp$k_opt, "\n")
  cat("SNP B:", res_snp$B, "\n")

  # -------------------------------------------------------
  # Grid
  # -------------------------------------------------------

  x1_grid <- seq(min(x1), max(x1), length.out = 50)
  x2_grid <- seq(min(x2), max(x2), length.out = 50)

  grid <- expand.grid(x1_grid, x2_grid)

  X1g <- matrix(grid[,1], 50, 50)
  X2g <- matrix(grid[,2], 50, 50)

  X_grid <- as.matrix(grid)

  # -------------------------------------------------------
  # DGCV surface
  # -------------------------------------------------------

  h_opt <- res_dgcv$h_opt_gcv

  W_grid <- construct_W(X, h_opt, X_grid)

  z_dgcv <- matrix(W_grid %*% y, 50, 50)

  # -------------------------------------------------------
  # SNP surface
  # -------------------------------------------------------

  h_start <- res_snp$h_start
  Yk      <- res_snp$y_k_minus_1_opt

  W_grid_snp <- construct_W(X, h_start, X_grid)

  z_snp <- matrix(W_grid_snp %*% Yk, 50, 50)

  # -------------------------------------------------------
  # Scale for plotting
  # -------------------------------------------------------

  y_plot <- y / 1000
  z_snp_plot  <- z_snp / 1000
  z_dgcv_plot <- z_dgcv / 1000

  # -------------------------------------------------------
  # Plot
  # -------------------------------------------------------

  par(mfrow = c(1,2))

  persp(
    x1_grid,
    x2_grid,
    z_snp_plot,
    theta = 40,
    phi = 30,
    col = "lightblue",
    shade = 0.4,
    ticktype = "detailed",
    xlab = "Median Income",
    ylab = "Housing Median Age",
    zlab = "House Value (×10^3)",
    main = "SNP Surface"
  )

  points(
    trans3d(x1, x2, y_plot,
      persp(
        x1_grid,
        x2_grid,
        z_snp_plot,
        plot = FALSE
      )
    ),
    pch = 16,
    cex = 0.3,
    col = rgb(0,0,0,0.2)
  )

  persp(
    x1_grid,
    x2_grid,
    z_dgcv_plot,
    theta = 40,
    phi = 30,
    col = "salmon",
    shade = 0.4,
    ticktype = "detailed",
    xlab = "Median Income",
    ylab = "Housing Median Age",
    zlab = "House Value (×10^3)",
    main = "DGCV Surface"
  )

  invisible(list(
    dgcv = res_dgcv,
    snp = res_snp,
    rmse_dgcv = rmse_dg,
    rmse_snp = rmse_sn
  ))
}
