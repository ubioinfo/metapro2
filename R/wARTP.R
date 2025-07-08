
#' Weighted ARTP with Adaptive Permutation
#'
#' Computes weighted ARTP p-value using adaptive two-stage permutation and estimates signal sparsity.
#'
#' @param p_values A numeric vector of p-values (0 to 1)
#' @param cor_mat Correlation matrix (optional). If NULL, assumes independence.
#' @param weight A numeric vector of weights (optional). Default = equal weight
#' @param B1 Number of permutations for stage 1 (default 2000)
#' @param B2 Number of permutations for stage 2 if needed (default 10000)
#' @param threshold p-value cutoff to trigger stage 2 (default 0.05)
#'
#' @return A list with:
#' \describe{
#'   \item{p_value}{Final adjusted p-value}
#'   \item{signal_sparsity}{Estimated proportion of influential signals}
#' }
#'
#' @examples
#' set.seed(1)
#' p <- runif(30) / 5
#' res <- wARTP(p)
#' res$p_value
#' res$signal_sparsity
#'
#' @export
wARTP <- function(p_values, cor_mat = NULL, weight = NULL, B1 = 2000, B2 = 10000, threshold = 0.05) {
  
  d <- length(p_values)
  
  if (!is.numeric(p_values) || any(p_values < 0 | p_values > 1)) {
    stop("p_values must be numeric and between 0 and 1")
  }
  if (is.null(weight)) {
    weight <- rep(1, d)
  } else if (length(weight) != d || !is.numeric(weight) || any(weight < 0)) {
    stop("weight must be a numeric vector of length equal to p_values with non-negative values")
  }
  if (is.null(cor_mat)) {
    cor_mat <- diag(d)
  } else if (!is.matrix(cor_mat) || nrow(cor_mat) != d || ncol(cor_mat) != d) {
    stop("cor_mat must be a square matrix with dimension equal to length of p_values")
  }
  if (!is.numeric(B1) || B1 < 1) {
    stop("B1 must be a positive integer")
  }
  if (!is.numeric(B2) || B2 < 1) {
    stop("B2 must be a positive integer")
  }
  if (!is.numeric(threshold) || threshold <= 0 || threshold >= 1) {
    stop("threshold must be between 0 and 1")
  }
  
  weights <- weight / mean(weight)
  
  res1 <- wARTP_cpp(p_values, cor_mat, weights, as.integer(B1))
  p1 <- res1$p_value
  ss <- res1$signal_sparsity
  
  if (p1 <= threshold) {
    res2 <- wARTP_cpp(p_values, cor_mat, weights, as.integer(B2))
    return(list(p_value = res2$p_value, signal_sparsity = res2$signal_sparsity))
  }
  
  return(list(p_value = p1, signal_sparsity = ss))
}
