#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <algorithm>
#include <vector>

// Main wARTP function: 정확한 qgamma만 사용
// [[Rcpp::export]]
Rcpp::List wARTP_cpp(const arma::vec& p_values, const arma::mat& cor_mat, const arma::vec& weights, int B) {
  int d = p_values.n_elem;
  
  arma::vec sorted_p = arma::sort(p_values);
  arma::vec sorted_w = weights(arma::sort_index(p_values));

  arma::vec qgamma_obs(d);
  for (int i = 0; i < d; ++i) {
    qgamma_obs(i) = R::qgamma(sorted_p(i), sorted_w(i), 2.0, false, false);
  }
  arma::vec stat_obs = arma::cumsum(qgamma_obs);

  arma::mat perm_stats(d, B);
  arma::mat chol_cor = arma::chol(cor_mat);

  for (int b = 0; b < B; ++b) {
    arma::vec z = chol_cor.t() * arma::randn<arma::vec>(d);
    arma::vec p_null = 2.0 * (1.0 - arma::normcdf(arma::abs(z)));

    arma::vec sorted_pn = arma::sort(p_null);
    arma::vec sorted_wn = weights(arma::sort_index(p_null));

    arma::vec qgamma_null(d);
    for (int i = 0; i < d; ++i) {
      qgamma_null(i) = R::qgamma(sorted_pn(i), sorted_wn(i), 2.0, false, false);
    }
    perm_stats.col(b) = arma::cumsum(qgamma_null);
  }

  arma::vec emp_pvals(d);
  for (int k = 0; k < d; ++k) {
    emp_pvals(k) = (arma::sum(stat_obs(k) <= perm_stats.row(k)) + 1.0) / (B + 1.0);
  }

  double min_pval = emp_pvals.min();
  int k_star = emp_pvals.index_min();
  double signal_sparsity = (k_star + 1.0) / d;

  arma::vec adaptive_null(B);
  for (int b = 0; b < B; ++b) {
    arma::vec temp_p(d);
    for (int k = 0; k < d; ++k) {
      temp_p(k) = (arma::sum(perm_stats(k, b) <= perm_stats.row(k)) + 1.0) / (B + 1.0);
    }
    adaptive_null(b) = temp_p.min();
  }

  double final_p = (arma::sum(min_pval >= adaptive_null) + 1.0) / (B + 1.0);

  return Rcpp::List::create(
    Rcpp::Named("p_value") = final_p,
    Rcpp::Named("signal_sparsity") = signal_sparsity
  );
}
