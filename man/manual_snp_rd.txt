\name{SNP}
\alias{SNP}
\title{Stepwise Noise Peeling for Nadaraya-Watson Regression}
\usage{
SNP(x, y, num_h_points = 50)
}
\arguments{
\item{x}{Numeric vector of predictor values (sorted).}
\item{y}{Numeric vector of response values corresponding to x.}
\item{num_h_points}{Integer, number of bandwidth candidates to evaluate in Phase I (default: 50).}
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
}
\examples{
# Generate sample data
set.seed(123)
n <- 100
x <- sort(runif(n, 0, 1))
y <- sin(2*pi*x) + rnorm(n, 0, 0.1)

# Apply SNP smoothing
result <- SNP(x, y)
plot(x, y, pch=16, col="gray")
lines(x, result$y_k_opt, col="red", lwd=2)
}