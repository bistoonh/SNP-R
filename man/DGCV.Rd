\name{DGCV}
\alias{DGCV}
\title{Direct Generalized Cross-Validation for Nadaraya-Watson Regression}
\usage{
DGCV(x, y, num_h_points = 50)
}
\arguments{
\item{x}{Numeric vector of predictor values (sorted).}
\item{y}{Numeric vector of response values corresponding to x.}
\item{num_h_points}{Integer, number of bandwidth candidates (default: 50).}
}
\value{
A list containing:
\describe{
  \item{y_h_opt}{Final smoothed output using optimal bandwidth}
  \item{h_opt_gcv}{Optimal bandwidth selected by GCV}
  \item{gcv_h}{GCV scores for all bandwidth candidates}
  \item{time_elapsed}{Execution time in seconds}
}
}
\description{
Implements Direct Generalized Cross-Validation (DGCV) for bandwidth selection
in Nadaraya-Watson regression. This is the reference method that SNP approximates.
}
\examples{
# Generate sample data
set.seed(123)
n <- 1000
x <- sort(runif(n, 0, 1))
y <- sin(2*pi*x) + rnorm(n, 0, 0.1)

# Apply DGCV
result <- DGCV(x, y)
plot(x, y, pch=16, col="gray")
lines(x, result$y_h_opt, col="blue", lwd=2)
}