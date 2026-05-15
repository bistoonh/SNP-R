# =========================================================
# metrics.R
# =========================================================

# ---------------------------------------------------------
# Root Mean Squared Error
# ---------------------------------------------------------

#' RMSE
#'
#' Compute Root Mean Squared Error.
#'
#' @param y_true True values
#' @param y_pred Predicted values
#'
#' @return Numeric scalar
#'
#' @export
rmse <- function(y_true, y_pred) {

  y_true <- as.numeric(y_true)
  y_pred <- as.numeric(y_pred)

  sqrt(mean((y_true - y_pred)^2))
}


# ---------------------------------------------------------
# Shifted Mean Absolute Percentage Error
# ---------------------------------------------------------

#' Shifted MAPE
#'
#' Compute shifted Mean Absolute Percentage Error.
#'
#' @param y_true True values
#' @param y_pred Predicted values
#'
#' @return Numeric scalar
#'
#' @export
mape_shift <- function(y_true, y_pred) {

  y_true <- as.numeric(y_true)
  y_pred <- as.numeric(y_pred)

  min_val <- min(min(y_true), min(y_pred))

  a <- max(0.0, 1.0 - min_val)

  100 * mean(
    abs((y_true - y_pred) / (y_true + a))
  )
}
