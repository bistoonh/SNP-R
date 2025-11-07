# SNP: Stepwise Noise Peeling for Nadaraya-Watson Regression

<!-- badges: start -->
[![R-CMD-check](https://github.com/yourusername/SNP/workflows/R-CMD-check/badge.svg)](https://github.com/bistoonh/SNP-R/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/SNP)](https://CRAN.R-project.org/package=SNP)
<!-- badges: end -->

The **SNP** package implements the Stepwise Noise Peeling algorithm that bypasses bandwidth selection in Nadaraya-Watson regression by using iterative smoothing. SNP provides a scalable alternative to Direct Generalized Cross-Validation (DGCV) by converting continuous bandwidth optimization into discrete iteration selection, dramatically reducing computational cost while maintaining statistical equivalence.

## Installation

You can install the development version of SNP from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("bistoonh/SNP-R", force = TRUE)
```

## Quick Start

```r
library(SNP)

# Generate sample data
set.seed(123)
n <- 2000
x <- sort(runif(n, 0, 1))
y <- sin(2*pi*x) + rnorm(n, 0, 0.1)

# Apply SNP smoothing with default parameters
snp_result <- SNP(x, y)

# Compare with traditional DGCV
dgcv_result <- DGCV(x, y)

# Plot results
plot(x, y, pch=16, col="gray", main="SNP vs DGCV Comparison")
lines(x, snp_result$y_k_opt, col="red", lwd=2)
lines(x, dgcv_result$y_h_opt, col="blue", lwd=2)
legend("topright", c("SNP", "DGCV"), col=c("red", "blue"), lwd=2)

# Performance comparison
cat("SNP time:", snp_result$time_elapsed, "seconds\n")
cat("DGCV time:", dgcv_result$time_elapsed, "seconds\n")
```

## Parameter Tuning

SNP provides two key parameters for balancing speed and accuracy:

```r
# Faster computation (fewer bandwidth candidates and slices)
snp_fast <- SNP(x, y, num_h_points = 25, num_slices = 25)

# More thorough search (more bandwidth candidates and slices)  
snp_thorough <- SNP(x, y, num_h_points = 50, num_slices = 100)

# Performance comparison
cat("Fast SNP time:", snp_fast$time_elapsed, "seconds\n")
cat("Thorough SNP time:", snp_thorough$time_elapsed, "seconds\n")
```

## Key Features

- **Fast**: Orders of magnitude faster than DGCV for large datasets
- **Accurate**: Statistically equivalent results to DGCV
- **Adaptive**: Automatically adjusts bandwidth through iterative process
- **Configurable**: Tunable parameters for speed vs accuracy trade-offs
- **Robust**: Handles edge cases and various data sizes
- **Well-documented**: Comprehensive help files and examples

## Algorithm Overview

SNP operates in two phases:

1. **Phase I**: Constructs a conservative initial bandwidth using random slices of data and lightweight GCV within each slice
   - `num_slices`: Controls number of random data slices (default: 60)
   - `num_h_points`: Controls bandwidth candidates per slice (default: 40)
   
2. **Phase II**: Fixes the smoothing operator and repeatedly applies it, selecting optimal iterations via discrete GCV

This reformulation preserves the adaptivity of GCV while converting costly continuous bandwidth search into lightweight discrete selection.

## Main Functions

- `SNP(x, y, num_h_points = 40, num_slices = 60)`: Main Stepwise Noise Peeling algorithm
- `DGCV(x, y, num_h_points = 50)`: Direct Generalized Cross-Validation (reference method)  
- `construct_W(x, h)`: Construct Gaussian kernel weight matrix

## Performance

For datasets with n > 1000, SNP typically shows:
- **Speed**: Orders of magnitude faster than DGCV
- **Accuracy**: < 1% difference in MSE compared to DGCV
- **Memory**: More efficient memory usage due to iterative approach
- **Scalability**: Parameter tuning allows adaptation to computational constraints

## Citation

If you use this package in your research, please cite:

```

```

## Contributing



## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Issues

Found a bug? Have a feature request? Please [open an issue](https://github.com/bistoonh/SNP-R/issues) on GitHub.
