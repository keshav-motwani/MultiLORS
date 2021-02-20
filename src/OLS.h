#ifndef OLSFIT
#define OLSFIT

#include <RcppArmadillo.h>

using namespace Rcpp;

arma::colvec OLS(const arma::mat & XtX, const arma::mat & X, const arma::colvec & Y);

#endif