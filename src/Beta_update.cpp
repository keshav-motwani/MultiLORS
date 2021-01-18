#include "Beta_update.h"
#include "objective.h"
#include "prox.h"

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

// [[Rcpp::export]]
arma::mat update_Beta(const List & Y_list, const List & X_list, const List & L_list, int q, const List & indices_list, const List & XtX_list, const List & XtY_list, const arma::mat & Beta_old, double lambda, double s_Beta, double s) {

  bool line_search = true;
  double shrinkage = 0.5;

  arma::mat Beta;
  arma::mat gradient = compute_gradient_Beta(X_list, L_list, q, indices_list, XtX_list, XtY_list, Beta_old);
  double g_old = evaluate_g(Y_list, X_list, L_list, indices_list, Beta_old);

  while(line_search & (s != s_Beta)) {

    Beta = l1_prox(Beta_old - (s * gradient), s * lambda);
    double g_new = evaluate_g(Y_list, X_list, L_list, indices_list, Beta);
    arma::mat difference = (Beta_old - Beta) / s;

    if (g_new > g_old - s * arma::accu(gradient % difference) + 0.5 * s * arma::accu(arma::pow(difference, 2))) {
      s = std::max(shrinkage * s, s_Beta);
    } else {
      line_search = false;
    }

  }

  return Beta;

}

// [[Rcpp::export]]
double compute_s_Beta(const List & XtX_list, int p, int q, const List & dataset_indices_list) {

  arma::vec eigenvalues(q);
  eigenvalues.zeros();

  for (R_xlen_t i = 0; i < q; i++) {

    arma::mat block(p, p);
    block.zeros();

    arma::uvec dataset_indices = dataset_indices_list[i];

    for (R_xlen_t j = 0; j < dataset_indices.n_elem; j++) {

      NumericMatrix XtX = XtX_list[j];

      arma::mat XtX_(XtX.begin(), XtX.nrow(), XtX.ncol(), false);

      block += XtX_;

    }

    arma::vec eigenvalues_block = arma::eig_sym(block);

    eigenvalues[i] = eigenvalues_block(eigenvalues_block.n_elem - 1);

  }

  return 1 / arma::max(eigenvalues);

}
