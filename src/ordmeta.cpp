
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List ordmeta_cpp(NumericVector pvals, int B = 20000) {
  int n = pvals.size();
  if (n == 0) return List::create(Named("p_value") = NA_REAL, Named("opt_ord") = NA_INTEGER);

  NumericVector ord_p = clone(pvals).sort();

  // Compute observed marginal p-values
  NumericVector margin_p(n);
  for (int r = 0; r < n; ++r) {
    margin_p[r] = R::pbeta(ord_p[r], r + 1, n - r, 1, 0);
  }
  double obs = Rcpp::min(margin_p);
  int opt_ord = Rcpp::which_min(margin_p) + 1; // R index starts at 1

  int count = 0;
  NumericVector u(n), beta(n);

  for (int b = 0; b < B; ++b) {
    for (int i = 0; i < n; ++i) u[i] = R::runif(0, 1);
    std::sort(u.begin(), u.end());

    for (int r = 0; r < n; ++r) {
      beta[r] = R::pbeta(u[r], r + 1, n - r, 1, 0);
    }
    if (obs > Rcpp::min(beta)) count++;
  }

  double pval = static_cast<double>(count) / B;
  return List::create(Named("p_value") = pval, Named("opt_ord") = opt_ord);
}
