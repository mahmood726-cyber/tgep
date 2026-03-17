\name{print.tgep}
\alias{print.tgep}
\title{Print Method for TGEP Objects}
\description{
  Prints a summary of a TGEP meta-analysis result, including the ensemble
  estimate, standard error, confidence interval, guard weights, and individual
  guard estimates.
}
\usage{
\method{print}{tgep}(x, \dots)
}
\arguments{
  \item{x}{An object of class \code{"tgep"} returned by \code{\link{tgep_meta}}.}
  \item{\dots}{Additional arguments (currently ignored).}
}
\value{
  Invisibly returns the input object \code{x}.
}
\seealso{
  \code{\link{tgep_meta}}, \code{\link{plot.tgep}}
}
\examples{
library(metafor)
yi <- c(0.2, 0.4, 0.1, 0.3, 0.5)
vi <- c(0.04, 0.06, 0.03, 0.05, 0.07)
res <- tgep_meta(yi, vi, n_boot = 0)
print(res)
}
