# =========================================================
# realdata_1d.R
# =========================================================

# ---------------------------------------------------------
# Load California housing data
# ---------------------------------------------------------

load_realdata_1d <- function() {

  data("housing", package = "snpregR")

  df <- na.omit(housing)

  x <- df$median_income

  y <- df$median_house_value

  idx <- order(x)

  list(
    x = x[idx],
    y = y[idx]
  )
}


# ---------------------------------------------------------
# Main experiment
# ---------------------------------------------------------

#' California Housing Example (1D)
#'
#' Fits SNP-NW and Direct-GCV NW using:
#' median_income -> median_house_value
#'
#' @param seed Random seed
#'
#' @return list
#'
#' @export
realdata_1d <- function(seed = 111) {

  dat <- load_realdata_1d()

  x <- dat$x

  y <- dat$y

  # -------------------------------------------------------
  # DGCV
  # -------------------------------------------------------

  set.seed(seed)

  res_dgcv <- nw_direct_gcv(
    x,
    y,
    num_h_points = 30,
    mode = "random"
  )

  yhat_dgcv <- as.numeric(
    res_dgcv$y_train_opt
  )

  # -------------------------------------------------------
  # SNP
  # -------------------------------------------------------

  set.seed(seed)

  res_snp <- nw_snp(
    x,
    y,
    num_h_points = 30,
    num_slices = 50
  )

  yhat_snp <- as.numeric(
    res_snp$y_k_opt
  )

  # -------------------------------------------------------
  # Metrics
  # -------------------------------------------------------

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
  # Plot
  # -------------------------------------------------------

  plot(
    x,
    y,
    pch = 16,
    cex = 0.4,
    col = rgb(0, 0, 0, 0.35),
    xlab = "Median Income",
    ylab = "Median House Value",
    main = "California Housing (1D)"
  )

  lines(
    x,
    yhat_snp,
    col = "#1f77b4",
    lwd = 2.6
  )

  lines(
    x,
    yhat_dgcv,
    col = "#d62728",
    lwd = 2.2
  )

  legend(
    "topleft",
    legend = c(
      "SNP",
      "DGCV"
    ),
    col = c(
      "#1f77b4",
      "#d62728"
    ),
    lwd = c(2.6, 2.2),
    bty = "n"
  )

  # -------------------------------------------------------
  # Return
  # -------------------------------------------------------

  invisible(list(
    dgcv = res_dgcv,
    snp = res_snp,
    rmse_dgcv = rmse_dg,
    rmse_snp = rmse_sn
  ))
}
