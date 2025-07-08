
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// wARTP_cpp
Rcpp::List wARTP_cpp(const arma::vec& p_values, const arma::mat& cor_mat, const arma::vec& weights, int B);
RcppExport SEXP _metapro2_wARTP_cpp(SEXP p_valuesSEXP, SEXP cor_matSEXP, SEXP weightsSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type p_values(p_valuesSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type cor_mat(cor_matSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(wARTP_cpp(p_values, cor_mat, weights, B));
    return rcpp_result_gen;
END_RCPP
}

// ordmeta_cpp
Rcpp::List ordmeta_cpp(Rcpp::NumericVector pvals, int B);
RcppExport SEXP _metapro2_ordmeta_cpp(SEXP pvalsSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::NumericVector pvals(pvalsSEXP);
    int B = Rcpp::as<int>(BSEXP);
    rcpp_result_gen = Rcpp::wrap(ordmeta_cpp(pvals, B));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_metapro2_wARTP_cpp", (DL_FUNC) &_metapro2_wARTP_cpp, 4},
    {"_metapro2_ordmeta_cpp", (DL_FUNC) &_metapro2_ordmeta_cpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_metapro2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
