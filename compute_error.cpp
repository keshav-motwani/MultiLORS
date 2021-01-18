#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends("RcppArmadillo")]]

// [[Rcpp::export]]
double compute_error(const List & X_list, const mat & Beta, const List & indices_list, const List & Y_list) {

  R_xlen_t K = X_list.size();

  double error = 0;

  for (R_xlen_t i = 0; i < K; i++) {

    NumericMatrix X = X_list[i];
    NumericMatrix Y = Y_list[i];
    uvec indices = indices_list[i];

    mat X_(X.begin(), X.nrow(), X.ncol(), false);
    mat Y_(Y.begin(), Y.nrow(), Y.ncol(), false);

    error = error + accu(pow(Y_ - X_ * Beta.cols(indices - 1), 2));

  }

  return error;
}

// [[Rcpp::export]]
arma::mat l1_prox(const arma::mat & matrix, double lambda, bool intercept) {

  arma::mat thresholded(matrix.n_rows, matrix.n_cols);
  thresholded.zeros();

  arma::uvec greater = arma::find(matrix > lambda);
  thresholded.elem(greater) = matrix.elem(greater) - lambda;

  arma::uvec less = arma::find(matrix < -lambda);
  thresholded.elem(less) = matrix.elem(less) + lambda;

  if (intercept) {

    thresholded.row(0) = matrix.row(0);

  }

  return thresholded;

}

// [[Rcpp::export]]
arma::mat nuclear_prox(const arma::mat & matrix, double gamma) {

  arma::mat U, V;
  arma::vec d;
  arma::svd_econ(U, d, V, matrix);

  arma::vec zero(d.size());
  zero.zeros();

  return U * arma::diagmat(arma::max(d - gamma, zero)) * V.t();

}


