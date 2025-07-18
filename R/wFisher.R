#' Weighted Fisher's Method for Combining P-values
#'
#' Computes a weighted version of Fisher's method using Gamma transformation.
#'
#' @param p A numeric vector of p-values (0 to 1)
#' @param weight A numeric vector of weights (optional). Default = equal weights
#' @param is.onetail Logical. If TRUE (default), assumes one-sided test; otherwise performs two-sided directional test using \code{eff.sign}
#' @param eff.sign A numeric vector of effect signs (+1 or -1), required if \code{is.onetail = FALSE}
#'
#' @return A list containing:
#' \describe{
#'   \item{p}{Combined p-value from weighted Fisher's method}
#'   \item{overall.eff.direction}{Direction of overall effect ("+" or "-") if two-sided test is used}
#' }
#'
#' @examples
#' set.seed(1)
#' p <- runif(10) / 5
#' wFisher(p)
#'
#' @export

# Weighted Fisher's method (one-tailed only)
wFisher <- function(p, weight = NULL) {
  if (is.null(weight)) {
    weight <- rep(1, length(p))
  }
  idx.na <- which(is.na(p))
  if (length(idx.na) > 0) {
    p <- p[-idx.na]
    weight <- weight[-idx.na]
  }
  NP <- length(p)
  NS <- length(weight)
  if (NP != NS) {
    stop("The length of p and weight vector must be identical.")
  }
  N <- NS
  Ntotal <- sum(weight)
  ratio <- weight / Ntotal
  Ns <- N * ratio
  
  G <- mapply(function(pi, wi) {
    qgamma(p = pi, shape = wi, scale = 2, lower.tail = FALSE)
  }, p, Ns)
  
  Gsum <- sum(G)
  resultP <- pgamma(q = Gsum, shape = N, scale = 2, lower.tail = FALSE)
  return(list(p = min(1, resultP)))
}