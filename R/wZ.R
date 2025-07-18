#' wZ function
#'
#' Description for wZ.
#'
#' @param ... Arguments passed to the function.
#'
#' @return Output of wZ.
#'
#' @examples
#' # Example usage of wZ()
#' wZ(...)
#'
#' @export


#' Weighted Z-test
#'
#' Combines p-values using the weighted Z method.
#'
#' @param p Vector of p-values
#' @param weight Vector of weights (e.g., sqrt(sample size))
#' @param is.onetail Whether to treat as one-tailed test
#' @param eff.sign Signs of effect sizes (for two-tailed)
#' @return Combined p-value and overall effect direction
#' @importFrom stats qnorm pnorm
#' @export
wZ <- function(p, weight, is.onetail = TRUE, eff.sign = NULL) {
  if (is.null(weight)) stop("Weight must be provided.")
  if (is.onetail) {
    z <- qnorm(1 - p)
  } else {
    z <- qnorm(p / 2, lower.tail = FALSE) * sign(eff.sign)
  }
  z_comb <- sum(weight * z)
  z_sd <- sqrt(sum(weight^2))
  z_final <- z_comb / z_sd
  p_comb <- 2 * (1 - pnorm(abs(z_final)))
  eff_dir <- sign(z_final)
  return(list(p = p_comb, overall.eff.direction = eff_dir))
}
