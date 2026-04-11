#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat rbf_kernel(const arma::vec& x1, const arma::vec& x2,
                     double length_scale, double variance) {
  int n1 = x1.n_elem;
  int n2 = x2.n_elem;
  arma::mat K(n1, n2);
  
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      double diff = x1[i] - x2[j];
      K(i, j) = variance * exp(-0.5 * (diff * diff) /
                                 (length_scale * length_scale));
    }
  }
  
  return K;
}
