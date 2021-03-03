#include "prox.h"

// [[Rcpp::export]]
arma::mat l1_prox(const arma::mat & matrix, double lambda) {

  arma::mat thresholded(matrix.n_rows, matrix.n_cols);
  thresholded.zeros();

  arma::uvec greater = arma::find(matrix > lambda);
  thresholded.elem(greater) = matrix.elem(greater) - lambda;

  arma::uvec less = arma::find(matrix < -lambda);
  thresholded.elem(less) = matrix.elem(less) + lambda;

  thresholded.row(0) = matrix.row(0);

  return thresholded;

}

// [[Rcpp::export]]
List nuclear_prox(const arma::mat & matrix, double gamma) {

  arma::mat U, V;
  arma::vec d;
  arma::svd_econ(U, d, V, matrix);

  arma::vec zero(d.size());
  zero.zeros();

  d = arma::max(d - gamma, zero);

  arma::mat L =  U * arma::diagmat(d) * V.t();

  double nuclear_norm_penalty = gamma * arma::accu(d);

  return List::create(Named("L") = wrap(L),
                      Named("nuclear_norm_penalty") = wrap(nuclear_norm_penalty));

}