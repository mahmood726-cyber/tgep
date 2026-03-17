\name{compare_tgep}
\alias{compare_tgep}
\title{Compare TGEP with Standard Meta-Analysis Estimators}
\description{
  Runs TGEP alongside REML and HKSJ (REML with Knapp-Hartung adjustment)
  estimators from \pkg{metafor} and returns a data frame comparing point
  estimates, standard errors, confidence intervals, and p-values.
}
\usage{
compare_tgep(yi, vi, n_boot = 200, temperature = 1.0, conf.level = 0.95)
}
\arguments{
  \item{yi}{Numeric vector of observed effect sizes.}
  \item{vi}{Numeric vector of corresponding sampling variances (must be
    strictly positive and the same length as \code{yi}).}
  \item{n_boot}{Number of bootstrap resamples for the TGEP standard error.
    Default is 200.}
  \item{temperature}{Softmax temperature for the TGEP ensemble.
    Default is 1.0.}
  \item{conf.level}{Confidence level for confidence intervals.
    Default is 0.95.}
}
\value{
  A \code{data.frame} with one row per method and columns:
  \item{Method}{Character string: \code{"TGEP"}, \code{"REML"}, or
    \code{"HKSJ"}.}
  \item{Estimate}{Pooled effect size estimate.}
  \item{SE}{Standard error.}
  \item{CI_Lower}{Lower confidence limit.}
  \item{CI_Upper}{Upper confidence limit.}
  \item{PValue}{Two-sided p-value.}

  If REML or HKSJ fitting fails, the corresponding row is omitted.
}
\examples{
library(metafor)

set.seed(42)
k <- 10
vi <- runif(k, 0.01, 0.10)
yi <- rnorm(k, 0.5, sqrt(vi))

tab <- compare_tgep(yi, vi)
print(tab)
}
\author{Mahmood Ul Hassan}
