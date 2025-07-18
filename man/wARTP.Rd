
\name{wARTP}
\alias{wARTP}
\title{Weighted ARTP Combination Test}
\description{Computes weighted ARTP p-value using Gamma-transformed statistics with adaptive permutation.}
\usage{
wARTP(p_values, cor_mat = NULL, weight = NULL, B1 = 2000, B2 = 10000, threshold = 0.05)
}
\arguments{
  \item{p_values}{A numeric vector of p-values (0 to 1)}
  \item{cor_mat}{Correlation matrix (optional). If NULL, assumes independence}
  \item{weight}{A numeric vector of weights (optional). Default = equal weights}
  \item{B1}{Number of permutations for stage 1}
  \item{B2}{Number of permutations for stage 2 when p-value < 0.05}
  \item{threshold}{p-value cutoff to trigger stage 2}
}
\value{
A list with:
\itemize{
  \item \code{p_value} Final p-value
  \item \code{signal_sparsity} Estimated proportion of influential signals
}
}
\examples{
set.seed(1)
p <- runif(20) / 2
res <- wARTP(p)
res\$p_value
res\$signal_sparsity
}
