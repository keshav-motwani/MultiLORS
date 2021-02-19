#include "OLS.h"

// [[Rcpp::export]]
arma::colvec OLS(const arma::mat & X, const arma::colvec & Y) {

  arma::colvec coef = arma::solve(X, Y);

  return coef;

}
