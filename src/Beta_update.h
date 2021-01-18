#ifndef BETA_UPDATE
#define BETA_UPDATE

#include <RcppArmadillo.h>

using namespace Rcpp;

arma::mat compute_gradient_Beta(const List & X_list, const List & L_list, int q, const List & indices_list, const List & XtX_list, const List & XtY_list, const arma::mat & Beta_old);
arma::mat update_Beta(const List & Y_list, const List & X_list, const List & L_list, int q, const List & indices_list, const List & XtX_list, const List & XtY_list, const arma::mat & Beta_old, double lambda, double s_Beta, double s);

#endif