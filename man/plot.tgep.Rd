\name{plot.tgep}
\alias{plot.tgep}
\title{Plot Method for TGEP Objects}
\description{
  Produces a two-panel diagnostic plot for a fitted \code{"tgep"} object.
  The left panel shows a bar chart of the three guard weights.  The right
  panel displays the individual guard estimates alongside the ensemble
  estimate with its confidence interval.
}
\usage{
\method{plot}{tgep}(x, \dots)
}
\arguments{
  \item{x}{An object of class \code{"tgep"} as returned by
    \code{\link{tgep_meta}}.}
  \item{\dots}{Additional arguments (currently unused).}
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

res <- tgep_meta(yi, vi)
plot(res)
}
\author{Mahmood Ul Hassan}
