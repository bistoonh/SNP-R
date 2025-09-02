#' Paper Examples for SNP Package
#' 
#' This file contains the examples used in the research paper demonstrating
#' the performance of SNP vs DGCV on different types of functions and real data.

library(SNP)

# Example 1: Stepwise Function -----------------------------------------------

#' Generate stepwise function with heteroscedastic noise
example_stepwise <- function() {
  
  # Define stepwise function
  stepwise_function <- function(x) {
    y <- numeric(length(x))
    
    # Define intervals (steps)
    idx1 <- which(x <= 20)
    idx2 <- which(x > 20 & x <= 35)
    idx3 <- which(x > 35 & x <= 45)
    idx4 <- which(x > 45)
    
    # Function values in each interval
    y[idx1] <- 2
    y[idx2] <- -1
    y[idx3] <- 3
    y[idx4] <- 0.5
    
    return(y)
  }
  
  # Generate data with heteroscedastic noise
  set.seed(2025)
  n <- 5000
  x <- runif(n, 0, 60)
  x <- sample(x)  # shuffle
  x <- sort(x)
  y_true <- stepwise_function(x)
  
  # Different noise variance in each interval
  noise_sd <- numeric(n)
  noise_sd[x <= 20] <- 0.2
  noise_sd[x > 20 & x <= 35] <- 0.8
  noise_sd[x > 35 & x <= 45] <- 0.1
  noise_sd[x > 45] <- 1.5
  
  noise <- rnorm(n, 0, noise_sd)
  y <- y_true + noise
  
  cat("=== Stepwise Function Example ===\n")
  cat("Dataset size:", n, "\n")
  
  # Apply SNP (suppress verbose output)
  snp_result <- suppressMessages(SNP(x, y))
  
  # Apply DGCV (suppress verbose output)
  dgcv_result <- suppressMessages(DGCV(x, y, num_h_points = 50))
  
  # Calculate RMSE against true function
  rmse_snp <- sqrt(mean((snp_result$y_k_opt - y_true)^2))
  rmse_dgcv <- sqrt(mean((dgcv_result$y_h_opt - y_true)^2))
  speedup <- dgcv_result$time_elapsed / snp_result$time_elapsed
  
  cat(sprintf("\n--- Results Summary ---\n"))
  cat(sprintf("SNP:  h_start=%.4f, k_opt=%d, RMSE=%.4f, Time=%.4fs\n", 
              snp_result$h_start, snp_result$k_opt, rmse_snp, snp_result$time_elapsed))
  cat(sprintf("DGCV: h_opt=%.4f, RMSE=%.4f, Time=%.4fs\n", 
              dgcv_result$h_opt_gcv, rmse_dgcv, dgcv_result$time_elapsed))
  cat(sprintf("Speedup: %.2fx\n", speedup))
  
  # Plot results (sample 2000 points for clarity)
  sample_idx <- sample(1:n, 2000)
  x_sample <- x[sample_idx]
  y_sample <- y[sample_idx]
  y_true_sample <- y_true[sample_idx]
  snp_sample <- snp_result$y_k_opt[sample_idx]
  dgcv_sample <- dgcv_result$y_h_opt[sample_idx]
  
  # Sort by x for proper line plotting
  sort_idx <- order(x_sample)
  x_sample <- x_sample[sort_idx]
  y_sample <- y_sample[sort_idx]
  y_true_sample <- y_true_sample[sort_idx]
  snp_sample <- snp_sample[sort_idx]
  dgcv_sample <- dgcv_sample[sort_idx]
  
  plot(x_sample, y_sample, pch=16, cex=0.3, col="gray70",
       main="Stepwise Function: SNP vs DGCV", 
       xlab="x", ylab="y", ylim=range(c(y_sample, y_true_sample)))
  lines(x_sample, y_true_sample, col="black", lwd=2, type="l")
  lines(x_sample, snp_sample, col="red", lwd=2)
  lines(x_sample, dgcv_sample, col="blue", lwd=2, lty=2)
  legend("topright", c("True Function", "SNP", "DGCV", "Data"), 
         col=c("black", "red", "blue", "gray70"), 
         lwd=c(2, 2, 2, NA), lty=c(1, 1, 2, NA), pch=c(NA, NA, NA, 16))
  
  return(list(
    performance = list(
      rmse_snp = rmse_snp,
      rmse_dgcv = rmse_dgcv
    )
  ))
}


# Example 2: Complex Wavy Function ------------------------------------------

#' Generate complex wavy function with mixed noise
example_wavy <- function() {
  
  set.seed(2025)
  n <- 10000
  x <- runif(n, 0, 40)
  x <- sample(x)
  x <- sort(x)
  
  # Complex function with multiple components
  f_true <- function(x) {
    y <- numeric(length(x))
    
    # Part 1: Sinusoidal with variable frequency and amplitude
    idx1 <- which(x <= 10)
    y[idx1] <- sin(0.3 * x[idx1]) * (1 + 0.3 * cos(0.5 * x[idx1]))
    
    # Part 2: Negative parabolic trend with steep slope
    idx2 <- which(x > 10 & x <= 20)
    y[idx2] <- -0.05 * (x[idx2] - 15)^2 + 0.5 * sin(0.7 * x[idx2])
    
    # Part 3: Linear with positive slope and severe noise
    idx3 <- which(x > 20 & x <= 30)
    y[idx3] <- 0.2 * (x[idx3] - 20) + 0.1 * sin(2 * x[idx3])
    
    # Part 4: High frequency oscillation with gradually decreasing amplitude
    idx4 <- which(x > 30)
    y[idx4] <- 0.5 * sin(5 * x[idx4]) * exp(-0.1 * (x[idx4] - 30))
    
    return(y)
  }
  
  y_true <- f_true(x)
  
  # Heterogeneous noise: mixed Gaussian and sparse noise
  noise_sd <- 0.2 + 0.8 * (x > 25) + 0.5 * sin(0.2 * x)^2
  noise_gauss <- rnorm(n, 0, noise_sd)
  noise_sparse <- rnorm(n, 0, 3) * (runif(n) < 0.02)  # 2% strong jumps
  noise <- noise_gauss + noise_sparse
  
  y <- y_true + noise
  
  cat("=== Complex Wavy Function Example ===\n")
  cat("Dataset size:", n, "\n")
  
  # Apply SNP (suppress verbose output)
  snp_result <- suppressMessages(SNP(x, y))
  
  # Apply DGCV (suppress verbose output)
  dgcv_result <- suppressMessages(DGCV(x, y, num_h_points = 50))
  
  # Calculate RMSE against true function
  rmse_snp <- sqrt(mean((snp_result$y_k_opt - y_true)^2))
  rmse_dgcv <- sqrt(mean((dgcv_result$y_h_opt - y_true)^2))
  speedup <- dgcv_result$time_elapsed / snp_result$time_elapsed
  
  cat(sprintf("\n--- Results Summary ---\n"))
  cat(sprintf("SNP:  h_start=%.4f, k_opt=%d, RMSE=%.4f, Time=%.4fs\n", 
              snp_result$h_start, snp_result$k_opt, rmse_snp, snp_result$time_elapsed))
  cat(sprintf("DGCV: h_opt=%.4f, RMSE=%.4f, Time=%.4fs\n", 
              dgcv_result$h_opt_gcv, rmse_dgcv, dgcv_result$time_elapsed))
  cat(sprintf("Speedup: %.2fx\n", speedup))
  
  # Plot results (sample 2000 points for clarity)
  sample_idx <- sample(1:n, 2000)
  x_sample <- x[sample_idx]
  y_sample <- y[sample_idx]
  y_true_sample <- y_true[sample_idx]
  snp_sample <- snp_result$y_k_opt[sample_idx]
  dgcv_sample <- dgcv_result$y_h_opt[sample_idx]
  
  # Sort by x for proper line plotting
  sort_idx <- order(x_sample)
  x_sample <- x_sample[sort_idx]
  y_sample <- y_sample[sort_idx]
  y_true_sample <- y_true_sample[sort_idx]
  snp_sample <- snp_sample[sort_idx]
  dgcv_sample <- dgcv_sample[sort_idx]
  
  plot(x_sample, y_sample, pch=16, cex=0.3, col="gray70",
       main="Complex Wavy Function: SNP vs DGCV", 
       xlab="x", ylab="y")
  lines(x_sample, y_true_sample, col="black", lwd=2)
  lines(x_sample, snp_sample, col="red", lwd=2)
  lines(x_sample, dgcv_sample, col="blue", lwd=2, lty=2)
  legend("topright", c("True Function", "SNP", "DGCV", "Data"), 
         col=c("black", "red", "blue", "gray70"), 
         lwd=c(2, 2, 2, NA), lty=c(1, 1, 2, NA), pch=c(NA, NA, NA, 16))
  
  return(list(
    performance = list(
      rmse_snp = rmse_snp,
      rmse_dgcv = rmse_dgcv
    )
  ))
}


# Example 3: Real Data - California Housing ---------------------------------

#' California Housing Dataset Example
example_california_housing <- function() {
  
  # Load California Housing data
  tryCatch({
    if (!requireNamespace("readr", quietly = TRUE)) {
      stop("Package 'readr' is required for this example. Please install it with: install.packages('readr')")
    }
    
    california_housing <- readr::read_csv(
      "https://raw.githubusercontent.com/ageron/handson-ml/master/datasets/housing/housing.csv",
      show_col_types = FALSE
    )
    
  }, error = function(e) {
    cat("Error loading California Housing dataset:\n")
    cat("Please check internet connection or install readr package\n")
    cat("Error message:", e$message, "\n")
    return(NULL)
  })
  
  if (exists("california_housing") && !is.null(california_housing)) {
    
    # Extract variables
    x <- california_housing$median_income
    y <- california_housing$median_house_value
    
    # Remove missing values
    complete_cases <- complete.cases(x, y)
    x <- x[complete_cases]
    y <- y[complete_cases]
    
    # Sort by x for better visualization
    sorted_idx <- order(x)
    x <- x[sorted_idx]
    y <- y[sorted_idx]
    
    n <- length(x)
    
    cat("=== California Housing Dataset Example ===\n")
    cat("Dataset size:", n, "\n")
    cat("Median income range: [", round(min(x), 2), ",", round(max(x), 2), "]\n")
    cat("House value range: [", round(min(y/1000), 0), "K,", round(max(y/1000), 0), "K]\n")
    
    # Apply SNP (suppress verbose output)
    snp_result <- suppressMessages(SNP(x, y))
    
    # Apply DGCV (suppress verbose output)
    dgcv_result <- suppressMessages(DGCV(x, y, num_h_points = 50))
    
    # Calculate RMSE (no true function, so use cross-method comparison)
    rmse_diff <- sqrt(mean((snp_result$y_k_opt - dgcv_result$y_h_opt)^2))
    speedup <- dgcv_result$time_elapsed / snp_result$time_elapsed
    
    cat(sprintf("\n--- Results Summary ---\n"))
    cat(sprintf("SNP:  h_start=%.4f, k_opt=%d, Time=%.4fs\n", 
                snp_result$h_start, snp_result$k_opt, snp_result$time_elapsed))
    cat(sprintf("DGCV: h_opt=%.4f, Time=%.4fs\n", 
                dgcv_result$h_opt_gcv, dgcv_result$time_elapsed))
    cat(sprintf("Speedup: %.2fx, RMSE difference: %.6f\n", speedup, rmse_diff))
    
    # Plot results (sample for better visualization)
    if (n > 3000) {
      sample_idx <- sample(1:n, 3000)
      x_plot <- x[sample_idx]
      y_plot <- y[sample_idx] / 1000  # Convert to thousands
      snp_plot <- snp_result$y_k_opt[sample_idx] / 1000
      dgcv_plot <- dgcv_result$y_h_opt[sample_idx] / 1000
    } else {
      x_plot <- x
      y_plot <- y / 1000
      snp_plot <- snp_result$y_k_opt / 1000
      dgcv_plot <- dgcv_result$y_h_opt / 1000
    }
    
    # Sort by x for proper line plotting
    sort_idx <- order(x_plot)
    x_plot <- x_plot[sort_idx]
    y_plot <- y_plot[sort_idx]
    snp_plot <- snp_plot[sort_idx]
    dgcv_plot <- dgcv_plot[sort_idx]
    
    plot(x_plot, y_plot, pch=16, cex=0.4, col="gray70",
         main="California Housing: Median Income vs House Value", 
         xlab="Median Income", ylab="Median House Value (thousands $)")
    lines(x_plot, snp_plot, col="red", lwd=2)
    lines(x_plot, dgcv_plot, col="blue", lwd=2, lty=2)
    legend("topright", c("SNP", "DGCV", "Data"), 
           col=c("red", "blue", "gray70"), 
           lwd=c(2, 2, NA), lty=c(1, 2, NA), pch=c(NA, NA, 16))
    
    return(list(
      performance = list(
        rmse_diff = rmse_diff
      )
    ))
    
  } else {
    cat("Could not load California Housing dataset. Skipping this example.\n")
    return(NULL)
  }
}


# Run All Examples Function -------------------------------------------------

#' Run all paper examples
run_all_examples <- function() {
  cat("Running all paper examples...\n\n")
  
  # Example 1: Stepwise
  result1 <- example_stepwise()
  readline(prompt = "Press [Enter] to continue to next example...")
  
  # Example 2: Wavy
  result2 <- example_wavy()
  readline(prompt = "Press [Enter] to continue to next example...")
  
  # Example 3: Real data
  result3 <- example_california_housing()
  
  cat("\n=== All Examples Complete ===\n")
  
  # Summary table
  if (!is.null(result1) && !is.null(result2) && !is.null(result3)) {
    cat("\nPerformance Summary:\n")
    cat(sprintf("%-15s %8s %8s %8s %8s\n", "Example", "SNP Time", "DGCV Time", "Speedup", "SNP RMSE"))
    cat(paste(rep("-", 50), collapse=""), "\n")
    cat(sprintf("%-15s %8.4f %8.4f %8.2fx %8.4f\n", 
                "Stepwise", result1$performance$snp_time, result1$performance$dgcv_time, 
                result1$performance$speedup, result1$performance$rmse_snp))
    cat(sprintf("%-15s %8.4f %8.4f %8.2fx %8.4f\n", 
                "Wavy", result2$performance$snp_time, result2$performance$dgcv_time, 
                result2$performance$speedup, result2$performance$rmse_snp))
    cat(sprintf("%-15s %8.4f %8.4f %8.2fx %8.6f\n", 
                "Housing", result3$performance$snp_time, result3$performance$dgcv_time, 
                result3$performance$speedup, result3$performance$rmse_diff))
  }
  
  return(list(
    stepwise = result1,
    wavy = result2,
    housing = result3
  ))
}