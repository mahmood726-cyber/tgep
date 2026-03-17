\name{plot_sparsity}
\alias{plot_sparsity}
\title{Plot Guard Weight Sparsity Across Temperatures}
\description{
  Evaluates the TGEP ensemble at multiple temperature values and plots how
  the three guard weights change as temperature varies.  Low temperatures
  concentrate weight on the best guard; high temperatures yield nearly
  equal weighting.
}
\usage{
plot_sparsity(yi, vi, temps = c(0.01, 0.1, 0.5, 1, 5, 10))
}
\arguments{
  \item{yi}{Numeric vector of observed effect sizes.}
  \item{vi}{Numeric vector of corresponding sampling variances (must be
    strictly positive and the same length as \code{yi}).}
  \item{temps}{Numeric vector of temperature values to evaluate.
    Default is \code{c(0.01, 0.1, 0.5, 1, 5, 10)}.}
}
\value{
  Called for its side effect of producing a plot.  Returns \code{NULL}
  invisibly.
}
\examples{
library(metafor)

set.seed(42)
k <- 8
vi <- runif(k, 0.01, 0.15)
yi <- rnorm(k, 0.3, sqrt(vi))

plot_sparsity(yi, vi)
}
\author{Mahmood Ul Hassan}
