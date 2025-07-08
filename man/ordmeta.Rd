
\name{ordmeta}
\alias{ordmeta}
\title{Ordinal Meta-analysis with Monte Carlo Estimation}
\description{Combines p-values using order-statistics and marginal Beta distributions with Monte Carlo approximation.}
\usage{
ordmeta(pvals, B = 20000)
}
\arguments{
  \item{pvals}{A numeric vector of p-values}
  \item{B}{Number of Monte Carlo replicates}
}
\value{
A list with:
\itemize{
  \item \code{p_value} Combined p-value based on order statistics
  \item \code{opt_ord} Optimal rank where the minimum marginal p-value occurs
}
}
\examples{
set.seed(123)
p <- runif(30) / 2
ordmeta(p)
}
