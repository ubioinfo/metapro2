// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// ordmeta_cpp
List ordmeta_cpp(NumericVector pvals, int B);
RcppExport SEXP _metapro2_ordmeta_cpp(SEXP pvalsSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pvals(pvalsSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(ordmeta_cpp(pvals, B));
    return rcpp_result_gen;
END_RCPP
}
// wARTP_cpp
Rcpp::List wARTP_cpp(const arma::vec& p_values, const arma::mat& cor_mat, Rcpp::Nullable<Rcpp::NumericVector> weights_opt, int B);
RcppExport SEXP _metapro2_wARTP_cpp(SEXP p_valuesSEXP, SEXP cor_matSEXP, SEXP weights_optSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type p_values(p_valuesSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type cor_mat(cor_matSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type weights_opt(weights_optSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(wARTP_cpp(p_values, cor_mat, weights_opt, B));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_metapro2_ordmeta_cpp", (DL_FUNC) &_metapro2_ordmeta_cpp, 2},
    {"_metapro2_wARTP_cpp", (DL_FUNC) &_metapro2_wARTP_cpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_metapro2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
