#include "OLS.h"

// [[Rcpp::export]]
arma::colvec OLS(const arma::mat & XtX, const arma::mat & X, const arma::colvec & Y) {

  arma::colvec coef = arma::solve(XtX, X.t() * Y);

  return coef;

}
