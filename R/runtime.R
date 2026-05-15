# =========================================================
# runtime.R
# =========================================================

# ---------------------------------------------------------
# Synthetic data generator
# ---------------------------------------------------------

generate_runtime_data <- function(
    n,
    d,
    noise_scale = 0.2
) {
  
  if (d == 1) {
    
    X1 <- runif(n)
    
    X <- matrix(X1, ncol = 1)
    
    y_true <- 1.5 * sin(2 * pi * X1)
    
  } else if (d == 2) {
    
    X1 <- runif(n)
    
    X2 <- runif(n, -3, 3)
    
    X <- cbind(X1, X2)
    
    y_true <-
      1.5 * sin(2 * pi * X1) +
      0.8 * cos(3 * X2)
    
  } else if (d == 3) {
    
    X1 <- runif(n)
    
    X2 <- runif(n, -3, 3)
    
    X3 <- runif(n, 0, 10)
    
    X <- cbind(X1, X2, X3)
    
    y_true <-
      1.5 * sin(2 * pi * X1) +
      0.8 * cos(3 * X2) +
      0.3 * (X3 - 5)^2 * exp(-0.2 * X3)
    
  } else {
    
    stop("d must be 1, 2, or 3")
    
  }
  
  noise <- rnorm(n, sd = noise_scale)
  
  y <- y_true + noise
  
  list(
    X = X,
    y = y,
    y_true = y_true
  )
}


# ---------------------------------------------------------
# Plot helper
# ---------------------------------------------------------

.make_runtime_plot <- function(
    df_mean
    
) {
  
  if (nrow(df_mean) == 0) {
    
    cat("[runtime] No results available.\n")
    
    return(invisible(NULL))
  }
  
  
  par(mfrow = c(2,2))
  
  methods <- unique(df_mean$method)
  
  dims_sorted <- sort(unique(df_mean$dim))
  
  # -------------------------------------------------------
  # Runtime
  # -------------------------------------------------------
  
  plot(
    NULL,
    xlim = range(df_mean$n),
    ylim = range(df_mean$time_elapsed),
    xlab = "n",
    ylab = "Time (seconds)",
    main = "Runtime vs Sample Size"
  )
  
  cols <- c("blue", "red", "darkgreen", "purple")
  
  idx_col <- 1
  
  legends <- c()
  
  for (d in dims_sorted) {
    
    for (method in methods) {
      
      df_sub <- df_mean[
        df_mean$dim == d &
          df_mean$method == method,
      ]
      
      df_sub <- df_sub[order(df_sub$n), ]
      
      lines(
        df_sub$n,
        df_sub$time_elapsed,
        type = "b",
        lwd = 2,
        col = cols[idx_col]
      )
      
      legends <- c(
        legends,
        paste0(method, " (d=", d, ")")
      )
      
      idx_col <- idx_col + 1
    }
  }
  
  legend(
    "topleft",
    legend = legends,
    col = cols[seq_along(legends)],
    lty = 1,
    pch = 1,
    cex = 0.8
  )
  
  grid()
  
  # -------------------------------------------------------
  # RMSE
  # -------------------------------------------------------
  
  plot(
    NULL,
    xlim = range(df_mean$n),
    ylim = range(df_mean$rmse),
    xlab = "n",
    ylab = "RMSE",
    main = "RMSE vs Sample Size"
  )
  
  idx_col <- 1
  
  for (d in dims_sorted) {
    
    for (method in methods) {
      
      df_sub <- df_mean[
        df_mean$dim == d &
          df_mean$method == method,
      ]
      
      df_sub <- df_sub[order(df_sub$n), ]
      
      lines(
        df_sub$n,
        df_sub$rmse,
        type = "b",
        lwd = 2,
        col = cols[idx_col]
      )
      
      idx_col <- idx_col + 1
    }
  }
  
  grid()
  
  # -------------------------------------------------------
  # MAPE
  # -------------------------------------------------------
  
  plot(
    NULL,
    xlim = range(df_mean$n),
    ylim = range(df_mean$mape_shift),
    xlab = "n",
    ylab = "MAPE (%)",
    main = "MAPE vs Sample Size"
  )
  
  idx_col <- 1
  
  for (d in dims_sorted) {
    
    for (method in methods) {
      
      df_sub <- df_mean[
        df_mean$dim == d &
          df_mean$method == method,
      ]
      
      df_sub <- df_sub[order(df_sub$n), ]
      
      lines(
        df_sub$n,
        df_sub$mape_shift,
        type = "b",
        lwd = 2,
        col = cols[idx_col]
      )
      
      idx_col <- idx_col + 1
    }
  }
  
  grid()
  
  dev.off()
  
  pdf(
    file.path(
      output_dir,
      "runtime_comparison.pdf"
    ),
    width = 11,
    height = 7
  )
  
  par(mfrow = c(2,2))
  
  plot(
    NULL,
    xlim = range(df_mean$n),
    ylim = range(df_mean$time_elapsed),
    xlab = "n",
    ylab = "Time (seconds)",
    main = "Runtime vs Sample Size"
  )
  
  idx_col <- 1
  
  for (d in dims_sorted) {
    
    for (method in methods) {
      
      df_sub <- df_mean[
        df_mean$dim == d &
          df_mean$method == method,
      ]
      
      df_sub <- df_sub[order(df_sub$n), ]
      
      lines(
        df_sub$n,
        df_sub$time_elapsed,
        type = "b",
        lwd = 2,
        col = cols[idx_col]
      )
      
      idx_col <- idx_col + 1
    }
  }
  
  legend(
    "topleft",
    legend = legends,
    col = cols[seq_along(legends)],
    lty = 1,
    pch = 1,
    cex = 0.8
  )
  
  grid()
  
  par(mfrow = c(1,1))
  
}


# ---------------------------------------------------------
# Runtime benchmark
# ---------------------------------------------------------

#' Runtime Scaling Benchmark
#'
#' @param n_list Sample sizes
#' @param dims Dimensions
#' @param n_rep Number of repetitions
#' @param seed Random seed
#'
#' @return list
#'
#' @export
runtime <- function(
    n_list = c(
      500,
      1500,
      3000,
      8000,
      13000,
      20000,
      30000
    ),
    dims = c(1, 2, 3),
    n_rep = 10,
    seed = 111
) {
  
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  records <- list()
  
  counter <- 1
  
  for (n in n_list) {
    
    cat(
      "\n========== Sample size n =",
      n,
      "==========\n"
    )
    
    for (d in dims) {
      
      cat(
        "[runtime] Starting dimension d =",
        d,
        "\n"
      )
      
      for (rep in seq_len(n_rep)) {
        
        current_seed <-
          seed +
          (n * 10000) +
          (d * 1000) +
          rep
        
        set.seed(current_seed)
        
        dat <- generate_runtime_data(
          n = n,
          d = d,
          noise_scale = 0.2
        )
        
        X <- dat$X
        y <- dat$y
        y_true <- dat$y_true
        
        cat(
          "  rep",
          rep,
          "/",
          n_rep,
          "(Seed:",
          current_seed,
          ")\n"
        )
        
        # ---------------------------------------------------
        # DGCV
        # ---------------------------------------------------
        
        out_dg <- nw_direct_gcv(
          X,
          y,
          num_h_points = 30,
          mode = "random"
        )
        
        yhat_dg <- as.numeric(
          out_dg$y_train_opt
        )
        
        records[[counter]] <- data.frame(
          dim = d,
          n = n,
          rep = rep,
          method = "DGCV",
          rmse = rmse(y_true, yhat_dg),
          mape_shift = mape_shift(
            y_true,
            yhat_dg
          ),
          time_elapsed = out_dg$time_elapsed
        )
        
        counter <- counter + 1
        
        # ---------------------------------------------------
        # SNP
        # ---------------------------------------------------
        
        out_snp <- nw_snp(
          X,
          y,
          num_h_points = 30,
          num_slices = 50
        )
        
        yhat_snp <- as.numeric(
          out_snp$y_k_opt
        )
        
        records[[counter]] <- data.frame(
          dim = d,
          n = n,
          rep = rep,
          method = "SNP",
          rmse = rmse(y_true, yhat_snp),
          mape_shift = mape_shift(
            y_true,
            yhat_snp
          ),
          time_elapsed = out_snp$time_elapsed
        )
        
        counter <- counter + 1
      }
      
      # -----------------------------------------------------
      # Save incremental progress
      # -----------------------------------------------------
      
      df_raw <- tmp$df_raw
      
      df_mean <- tmp$df_mean
    }
    
  }
  return(list(df_raw = df_raw, df_mean = df_mean))
}
