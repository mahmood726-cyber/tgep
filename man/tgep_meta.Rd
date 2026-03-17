\name{tgep_meta}
\alias{tgep_meta}
\title{Triple-Guard Ensemble Pooling for Meta-Analysis}
\description{
  Computes a pooled effect size estimate using three specialized guard

  estimators (GRMA, WRD, SWA) integrated via leave-one-out cross-validation
  stacking with softmax weighting controlled by a temperature parameter.
  Bootstrap resampling provides standard errors that account for model
  uncertainty.
}
\usage{
tgep_meta(yi, vi, n_boot = 200, temperature = 1.0, n_cores = 1,
          conf.level = 0.95)
}
\arguments{
  \item{yi}{Numeric vector of observed effect sizes (e.g., log odds ratios,
    standardized mean differences).}
  \item{vi}{Numeric vector of corresponding sampling variances (must be
    strictly positive and the same length as \code{yi}).}
  \item{n_boot}{Number of bootstrap resamples for standard error estimation.
    Set to 0 to use an analytical approximation instead.  Default is 200.}
  \item{temperature}{Softmax temperature controlling the sparsity of the
    ensemble weights.  Lower values concentrate weight on the best-performing
    guard; higher values spread weight more evenly.  Default is 1.0.}
  \item{n_cores}{Number of CPU cores for parallel bootstrapping (via
    \code{\link[parallel]{mclapply}}).  Parallelism is only used on non-Windows
    platforms; on Windows this argument is ignored.  Default is 1.}
  \item{conf.level}{Confidence level for the confidence interval.
    Default is 0.95.}
}
\value{
  An object of class \code{"tgep"}, which is a list with components:
  \item{estimate}{Pooled effect size estimate.}
  \item{se}{Standard error of the estimate (bootstrap or analytical).}
  \item{ci_lb}{Lower bound of the confidence interval.}
  \item{ci_ub}{Upper bound of the confidence interval.}
  \item{pvalue}{Two-sided p-value based on the Wald test.}
  \item{weights}{Named numeric vector of ensemble weights for the three
    guards (GRMA, WRD, SWA).}
  \item{guards}{Named numeric vector of individual guard point estimates.}
  \item{k}{Number of studies.}
  \item{temp}{Temperature value used.}
  \item{conf.level}{Confidence level used.}
}
\examples{
library(metafor)

# Simulate a small meta-analysis
set.seed(42)
k <- 8
theta <- 0.3
vi <- runif(k, 0.01, 0.15)
yi <- rnorm(k, theta, sqrt(vi))

# Fit TGEP
res <- tgep_meta(yi, vi)
print(res$estimate)
print(res$weights)

# Without bootstrap (analytical SE)
res0 <- tgep_meta(yi, vi, n_boot = 0)
}
\author{Mahmood Ul Hassan}
