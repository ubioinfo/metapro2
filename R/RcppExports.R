# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

ordmeta_cpp <- function(pvals, B = 20000L) {
    .Call(`_metapro2_ordmeta_cpp`, pvals, B)
}

wARTP_cpp <- function(p_values, cor_mat, weights_opt, B) {
    .Call(`_metapro2_wARTP_cpp`, p_values, cor_mat, weights_opt, B)
}

