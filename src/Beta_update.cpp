#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat compute_gradient_Beta(const List & X_list, const List & L_list, int q, const List & indices_list, const List & XtX_list, const List & XtY_list, const arma::mat & Beta_old) {

  R_xlen_t K = X_list.size();

  int p = Beta_old.n_rows;

  arma::mat gradient(p, q);
  gradient.zeros();

  for (R_xlen_t i = 0; i < K; i++) {

    NumericMatrix XtX = XtX_list[i];
    NumericMatrix L = L_list[i];
    NumericMatrix XtY = XtY_list[i];
    NumericMatrix X = X_list[i];
    arma::uvec indices = indices_list[i];

    arma::mat XtX_(XtX.begin(), XtX.nrow(), XtX.ncol(), false);
    arma::mat L_(L.begin(), L.nrow(), L.ncol(), false);
    arma::mat XtY_(XtY.begin(), XtY.nrow(), XtY.ncol(), false);
    arma::mat X_(X.begin(), X.nrow(), X.ncol(), false);

    gradient.cols(indices - 1) += XtX_ * Beta_old.cols(indices - 1) + X_.t() * L_ - XtY_;

  }

  return gradient;
}
