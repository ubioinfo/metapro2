#' lancaster function
#'
#' Description for lancaster.
#'
#' @param ... Arguments passed to the function.
#'
#' @return Output of lancaster.
#'
#' @examples
#' # Example usage of lancaster()
#' lancaster(...)
#'
#' @export


#' Lancaster's method
#'
#' Combines p-values using Lancaster's weighted method.
#'
#' @param p Vector of p-values
#' @param df Degrees of freedom for each test
#' @return Combined p-value
#' @importFrom stats qchisq pchisq
#' @export
lancaster <- function(p, df) {
  chisq_values <- qchisq(1 - p, df = df)
  total_stat <- sum(chisq_values)
  p_combined <- pchisq(total_stat, df = sum(df), lower.tail = FALSE)
  return(p_combined)
}
