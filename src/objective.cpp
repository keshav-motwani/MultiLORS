#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
double evaluate_g(const List & Y_list, const List & X_list, const List & L_list, const List & indices_list, const arma::mat & Beta) {

  R_xlen_t K = X_list.size();

  double error = 0;

  for (R_xlen_t i = 0; i < K; i++) {

    NumericMatrix X = X_list[i];
    NumericMatrix Y = Y_list[i];
    NumericMatrix L = L_list[i];
    arma::uvec indices = indices_list[i];

    arma::mat X_(X.begin(), X.nrow(), X.ncol(), false);
    arma::mat Y_(Y.begin(), Y.nrow(), Y.ncol(), false);
    arma::mat L_(L.begin(), L.nrow(), L.ncol(), false);

    error += arma::accu(arma::pow(Y_ - X_ * Beta.cols(indices - 1) - L_, 2));

  }

  return error / 2;

}

// [[Rcpp::export]]
double l1_penalty(const arma::mat & Beta, double lambda) {

  int p = Beta.n_rows;

  return lambda * arma::accu(arma::abs(Beta.rows(1, p - 1)));

}

// [[Rcpp::export]]
double nuclear_norm_penalty(const List & L_list, double gamma, const arma::vec & gamma_weights) {

  R_xlen_t K = L_list.size();

  double penalty = 0;

  for (R_xlen_t i = 0; i < K; i++) {

    NumericMatrix L = L_list[i];

    arma::mat L_(L.begin(), L.nrow(), L.ncol(), false);

    arma::mat U, V;
    arma::vec d;
    arma::svd_econ(U, d, V, L_, "right");

    penalty += arma::accu(d) * gamma * gamma_weights[i];

  }

  return penalty;

}

// [[Rcpp::export]]
double evaluate_objective(const List & Y_list, const List & X_list, const List & L_list, const List & indices_list, const arma::mat & Beta, double lambda, double gamma, const arma::vec & gamma_weights) {

  double g = evaluate_g(Y_list, X_list, L_list, indices_list, Beta);
  double l1 = l1_penalty(Beta, lambda);
  double nuclear = nuclear_norm_penalty(L_list, gamma, gamma_weights);

  return g + l1 + nuclear;

}