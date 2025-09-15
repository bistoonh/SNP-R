\name{SNP}
\alias{SNP}
\title{Stepwise Noise Peeling for Nadaraya-Watson Regression}
\usage{
SNP(x, y, num_h_points = 40, num_slices = 60)
}
\arguments{
\item{x}{Numeric vector of predictor values (sorted).}

\item{y}{Numeric vector of response values corresponding to x.}

\item{num_h_points}{Integer, number of bandwidth candidates to evaluate
  within each slice (default: 40).}

\item{num_slices}{Integer, number of random slices to use for initial
  bandwidth estimation (default: 60).}
}
\value{
A list containing:
\describe{
  \item{y_k_opt}{Final smoothed output vector}
  \item{h_start}{Final chosen initial bandwidth}
  \item{k_opt}{Optimal number of iterations}
  \item{gcv_approx_k}{GCV values for each iteration}
  \item{time_elapsed}{Execution time in seconds}
}
}
\description{
Implements the Stepwise Noise Peeling (SNP) algorithm for efficient bandwidth
selection in Nadaraya-Watson regression with Gaussian kernels.
}
\details{
The SNP algorithm operates in two phases:
Phase I constructs a conservative initial bandwidth using random slices of data.
Phase II fixes the smoothing operator and repeatedly applies it, selecting optimal iterations via discrete GCV.

The \code{num_h_points} parameter controls the number of bandwidth candidates
tested within each slice during Phase I. Higher values provide more thorough
optimization but increase computational cost.

The \code{num_slices} parameter determines how many random data slices are used
for initial bandwidth estimation. More slices generally improve robustness of
the median aggregation but require additional computation.
}
\examples{
# Generate sample data
set.seed(123)
n <- 1500
x <- sort(runif(n, 0, 1))
y <- sin(2*pi*x) + rnorm(n, 0, 0.1)

# Apply SNP smoothing with default parameters
result <- SNP(x, y)
plot(x, y, pch=16, col="gray")
lines(x, result$y_k_opt, col="red", lwd=2)

# Apply SNP with custom parameters for faster computation
result_fast <- SNP(x, y, num_h_points = 20, num_slices = 30)
lines(x, result_fast$y_k_opt, col="blue", lwd=2, lty=2)

# Apply SNP with more thorough parameter search
result_thorough <- SNP(x, y, num_h_points = 60, num_slices = 100)
lines(x, result_thorough$y_k_opt, col="green", lwd=2, lty=3)

legend("topright", 
       c("Default", "Fast", "Thorough"), 
       col=c("red", "blue", "green"), 
       lwd=2, lty=c(1, 2, 3))
}