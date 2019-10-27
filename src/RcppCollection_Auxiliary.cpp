#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/*
 * 1. image normalization
 */
// [[Rcpp::export]]
arma::mat rcpp_01normalize(arma::mat input){
  const int m = input.n_rows;
  const int n = input.n_cols;

  arma::mat output(m,n,fill::zeros);
  const double a = input.min();
  const double b = input.max();
  const double denominator = 1/(b-a);
  for (int i=0;i<m;i++){
    for (int j=0;j<n;j++){
      output(i,j) = (input(i,j)-a)*denominator;
    }
  }
  return(output);
}
