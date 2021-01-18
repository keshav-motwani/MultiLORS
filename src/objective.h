#ifndef OBJECTIVE
#define OBJECTIVE

#include <RcppArmadillo.h>

using namespace Rcpp;

double evaluate_g(const List & Y_list, const List & X_list, const List & L_list, const List & indices_list, const arma::mat & Beta);
double l1_penalty(const arma::mat & Beta, double lambda);
double nuclear_norm_penalty(const List & L_list, double gamma, const arma::vec & gamma_weights);
double evaluate_objective(const List & Y_list, const List & X_list, const List & L_list, const List & indices_list, const arma::mat & Beta, double lambda, double gamma, const arma::vec & gamma_weights);

#endif
