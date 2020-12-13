#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//[[Rcpp::depends("RcppArmadillo")]]
//[[Rcpp::export]]
List svd_c(arma::mat X){

    arma::mat U;
    arma::vec s;
    arma::mat V;
    svd(U, s, V, X);

    List ret;


    return(ret);

}