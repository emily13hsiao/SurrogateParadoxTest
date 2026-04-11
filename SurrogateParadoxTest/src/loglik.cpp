#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double loglik_contrib(const arma::mat& Ck,
                      const arma::vec& yk,
                      const arma::vec& mk) {
  arma::vec v = yk - mk;
  arma::vec sol = arma::solve(Ck, v);
  double quad = arma::as_scalar(v.t() * sol);
  
  double logdet_val = 0.0;
  double sign = 0.0;
  arma::log_det(logdet_val, sign, Ck);
  
  return 0.5 * (logdet_val + quad);
}
