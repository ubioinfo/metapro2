#include <RcppArmadillo.h>
#include <algorithm> // For std::lower_bound and std::distance

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
Rcpp::List wARTP_cpp(const arma::vec& p_values,
                               const arma::mat& cor_mat,
                               Rcpp::Nullable<Rcpp::NumericVector> weights_opt,
                               int B) {
  int d = p_values.n_elem;
  arma::uvec sort_idx = arma::sort_index(p_values);
  arma::vec sorted_p = p_values(sort_idx);

  bool use_weights = false;
  arma::vec weights(d, arma::fill::ones);
  if (weights_opt.isNotNull()) {
    Rcpp::NumericVector w_input(weights_opt);
    if (w_input.size() == d) {
      weights = Rcpp::as<arma::vec>(w_input);
      if (arma::stddev(weights) > 1e-8) {
        use_weights = true;
      }
    }
  }

  arma::mat chol_cor = arma::chol(cor_mat);
  arma::mat perm_stats(d, B);
  arma::vec stat_obs(d);

  // This section remains the same: generate observed and permuted statistics
  if (use_weights) {
    arma::vec sorted_w = weights(sort_idx);
    for (int i = 0; i < d; ++i) {
      stat_obs(i) = R::qgamma(sorted_p(i), sorted_w(i), 2.0, false, false);
    }
    stat_obs = arma::cumsum(stat_obs);

    for (int b = 0; b < B; ++b) {
      arma::vec z = chol_cor.t() * arma::randn(d);
      arma::vec p_null = 2.0 * (1.0 - arma::normcdf(arma::abs(z)));
      arma::uvec sort_idxn = arma::sort_index(p_null);
      arma::vec sorted_pn = p_null(sort_idxn);
      arma::vec sorted_wn = weights(sort_idxn);
      arma::vec stat_null(d);
      for (int i = 0; i < d; ++i) {
        stat_null(i) = R::qgamma(sorted_pn(i), sorted_wn(i), 2.0, false, false);
      }
      perm_stats.col(b) = arma::cumsum(stat_null);
    }
  } else {
    for (int i = 0; i < d; ++i) {
      stat_obs(i) = -2.0 * std::log(sorted_p(i));
    }
    stat_obs = arma::cumsum(stat_obs);

    for (int b = 0; b < B; ++b) {
      arma::vec z = chol_cor.t() * arma::randn(d);
      arma::vec p_null = 2.0 * (1.0 - arma::normcdf(arma::abs(z)));
      arma::vec sorted_pn = arma::sort(p_null);
      arma::vec stat_null(d);
      for (int i = 0; i < d; ++i) {
        stat_null(i) = -2.0 * std::log(sorted_pn(i));
      }
      perm_stats.col(b) = arma::cumsum(stat_null);
    }
  }

  // Empirical p-values calculation (O(d*B), already reasonably efficient)
  arma::vec emp_pvals(d);
  for (int k = 0; k < d; ++k) {
    emp_pvals(k) = (arma::sum(stat_obs(k) <= perm_stats.row(k)) + 1.0) / (B + 1.0);
  }

  double min_pval = emp_pvals.min();
  int k_star = emp_pvals.index_min();
  double signal_sparsity = (k_star + 1.0) / d;

  // --- OPTIMIZED SECTION ---
  // Adaptive null distribution calculation
  // The original O(B^2 * d) nested loop was the main bottleneck.
  // The new approach reduces complexity to O(d * B * log(B)) while producing
  // the exact same result, leading to a significant speed-up.

  arma::vec adaptive_null(B);
  arma::mat p_val_matrix(d, B);

  // 1. For each rank (row k), calculate the p-value of each permuted statistic
  //    against its own null distribution (the other stats in the same row).
  for (int k = 0; k < d; ++k) {
    arma::rowvec row_k = perm_stats.row(k);
    arma::rowvec sorted_row_k = arma::sort(row_k); // Sort the row: O(B log B)

    for (int b = 0; b < B; ++b) {
      double val = row_k(b);
      // Use binary search (std::lower_bound) to efficiently find the number of
      // permuted stats greater than or equal to the current value 'val'.
      // This is much faster than a linear scan. Complexity: O(log B)
      auto it = std::lower_bound(sorted_row_k.begin(), sorted_row_k.end(), val);
      double count = std::distance(it, sorted_row_k.end());
      p_val_matrix(k, b) = (count + 1.0) / (B + 1.0);
    }
  }

  // 2. Find the minimum p-value across all ranks for each permutation.
  for (int b = 0; b < B; ++b) {
    adaptive_null(b) = p_val_matrix.col(b).min();
  }
  // --- END OF OPTIMIZED SECTION ---

  // Final p-value calculation
  double final_p = (arma::sum(min_pval >= adaptive_null) + 1.0) / (B + 1.0);

  return Rcpp::List::create(
    Rcpp::Named("p_value") = final_p,
    Rcpp::Named("signal_sparsity") = signal_sparsity,
    Rcpp::Named("used_weights") = use_weights
  );
}
