#ifndef PROX
#define PROX

#include <RcppArmadillo.h>

using namespace Rcpp;

arma::mat l1_prox(const arma::mat & matrix, double lambda);
List nuclear_prox(const arma::mat & matrix, double gamma);

#endif