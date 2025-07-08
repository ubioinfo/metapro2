
#' Ordinal Meta-analysis with Monte Carlo Estimation
#'
#' Computes order-statistics based combined p-value using marginal Beta distributions and Monte Carlo approximation.
#'
#' @param pvals A numeric vector of p-values (0 to 1)
#' @param B Number of Monte Carlo replicates (default 20000)
#'
#' @return A list with:
#' \describe{
#'   \item{p_value}{Combined p-value based on order statistics}
#'   \item{opt_ord}{Optimal rank (position) where the minimum marginal p-value occurs}
#' }
#'
#' @examples
#' set.seed(1)
#' p <- runif(25)
#' res <- ordmeta(p)
#' res$p_value
#' res$opt_ord
#'
#' @export
ordmeta <- function(pvals, B = 20000) {
  res <- ordmeta_cpp(pvals, B)
  return(list(p_value = res$p_value, opt_ord = res$opt_ord))
}
