# snpreg: Stepwise Noise Peeling for Nonparametric Regression

An R implementation of the Stepwise Noise Peeling (SNP) method for efficient bandwidth selection in Nadaraya-Watson kernel regression.

## Overview

SNP is a computationally efficient alternative to Generalized Cross-Validation (GCV) for bandwidth selection in nonparametric regression. It transforms the continuous $d$-dimensional bandwidth optimization problem into a discrete one-dimensional path, achieving significant speedups (over $30\times$) while maintaining or improving prediction accuracy.

### Key Features

- **Fast bandwidth selection**: Over $30\times$ faster than traditional GCV
- **High accuracy**: Comparable or better prediction performance than GCV
- **Scalability**: Efficient for high-dimensional problems
- **Stability**: Avoids local minima in bandwidth optimization
- **Simple API**: Easy-to-use interface for both 1D and multidimensional regression

### How It Works

SNP operates in two phases:

**Phase I (Coarse Search):**
- Evaluates bandwidths on a logarithmic grid using data slices
- Uses fast spectral approximation for effective degrees of freedom
- Identifies promising bandwidth regions efficiently

**Phase II (Refinement):**
- Performs adaptive restarts from Phase I candidates
- Uses approximate GCV on full data for fine-tuning
- Selects optimal bandwidth with minimal computational cost

## Installation

Currently, install directly from GitHub:
```r
# Install devtools if needed
if (!require("devtools")) install.packages("devtools")

# Install snpreg
devtools::install_github("bistoonh/SNP-R", force = TRUE)
```
### Requirements

- R >= 3.5.0
- stats (base R)
- graphics (base R)

## Quick Start

This example demonstrates basic usage of SNP for 1D regression:

```r
library(snpreg)

# Generate synthetic data
set.seed(111)
n <- 10000
X <- runif(n, 0, 10)
Y <- sin(2*X) + rnorm(n, 0, 0.25)

# Fit using SNP (fast)
result <- nw_snp(X, Y)
hstart <- result$h_start
kopt <- result$k_opt
B <- result$B
timeElapsed <- result$time_elapsed

cat(sprintf("h_start: %s, k_opt: %s, B: %s, time_elapsed: %.2f(Sec)\n", 
hstart, kopt, B, timeElapsed))
```

## Examples

### Example 1: 1D Smoothing with Comparison

This example compares SNP with Direct GCV on a noisy sine wave, demonstrating SNP's speed advantage while maintaining accuracy:

```r
library(snpreg)

# Generate noisy sine wave
set.seed(111)
n <- 2000

X <- sort(runif(n, 0, 4*pi))
Y_true <- sin(2*X)
Y <- Y_true + rnorm(n, 0, 0.3)

# Fit with SNP
result_snp <- nw_snp(X, Y, num_h_points = 30, num_slices = 50)

# Fit with Direct GCV (for comparison)
result_gcv <- nw_direct_gcv(X, Y, num_h_points = 30)

# Evaluate
rmse_snp <- rmse(Y_true, result_snp$y_k_opt)
rmse_gcv <- rmse(Y_true, result_gcv$y_train_opt)
mape_snp <- mape_shift(Y_true, result_snp$y_k_opt)
mape_gcv <- mape_shift(Y_true, result_gcv$y_train_opt)

cat(sprintf("SNP: Elapsed Time: %.1f(Sec), RMSE: %.3f, MAPE: %.2f\n", 
            result_snp$time_elapsed, rmse_snp, mape_snp))
cat(sprintf("GCV: Elapsed Time: %.1f(Sec), RMSE: %.3f, MAPE: %.2f\n", 
            result_gcv$time_elapsed, rmse_gcv, mape_gcv))

# Plot results
par(mfrow = c(1, 2))

plot(X, Y, pch = 16, cex = 0.5, col = 'gray', 
     xlab = "X", ylab = "y", 
     main = sprintf("SNP: RMSE=%.3f, MAPE=%.2f", rmse_snp, mape_snp))
lines(X, Y_true, col = "black", lwd = 2)
lines(X, result_snp$y_k_opt, col = "red", lwd = 2)
legend("topright", c("True function", "Training data", "SNP"), 
       col = c("black", "gray", "red"), lty = c(1, NA, 1), 
       pch = c(NA, 16, NA), lwd = c(2, NA, 2))
grid()

plot(X, Y, pch = 16, cex = 0.5, col = 'gray', 
     xlab = "X", ylab = "y", 
     main = sprintf("Direct GCV: RMSE=%.3f, MAPE=%.2f", rmse_gcv, mape_gcv))
lines(X, Y_true, col = "black", lwd = 2)
lines(X, result_gcv$y_train_opt, col = "green", lwd = 2)
legend("topright", c("True function", "Training data", "GCV"), 
       col = c("black", "gray", "green"), lty = c(1, NA, 1), 
       pch = c(NA, 16, NA), lwd = c(2, NA, 2))
grid()

par(mfrow = c(1, 1))
```

## Experiments Section

The package includes reproducible experiments from the paper. These functions run comprehensive benchmarks and generate publication-quality figures:

### Run runtime benchmark
```r
library(snpreg)

res <- runtime(n_list = c(500), dims = c(1, 2, 3))
res$mean_results```
### Run mixture experiment (synthetic data)
```r
library(snpreg)

mixture_experiment()
```
### Run 1D real data experiment
```r
library(snpreg)

realdata_1d()
```
### Run 2D real data experiment
```r
library(snpreg)

realdata_2d()
```
## Parameter Tuning

SNP provides two key parameters for balancing speed and accuracy:

```r
library(snpreg)

# Generate noisy sine wave
set.seed(111)
n <- 13000

X <- runif(n, 0, 4*pi)
Y_true <- sin(2*X)
Y <- Y_true + rnorm(n, 0, 0.3)

# Faster computation (fewer bandwidth candidates and slices)
result_fast <- nw_snp(X, Y, num_h_points = 20, num_slices = 25)

# More thorough search (more bandwidth candidates and slices)
result_thorough <- nw_snp(X, Y, num_h_points = 50, num_slices = 100)

# Performance comparison
rmse_fast <- rmse(Y_true, result_fast$y_k_opt)
rmse_thorough <- rmse(Y_true, result_thorough$y_k_opt)

cat(sprintf("Fast SNP time: %.2f seconds, RMSE: %.3f\n", 
result_fast$time_elapsed, rmse_fast))
cat(sprintf("Thorough SNP time: %.2f seconds, RMSE: %.3f\n", 
result_thorough$time_elapsed, rmse_thorough))
```

## License

MIT License - see LICENSE file for details.

## Author

**Bistoon Hosseini**  
Email: bistoon.hosseini@gmail.com  
GitHub: [github.com/bistoonh/SNP-R](https://github.com/bistoonh/SNP-R)

