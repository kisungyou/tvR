#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// auxiliary functions for 1d signal
arma::rowvec rcpp_diff(arma::rowvec x){
  const int n = x.n_elem;
  arma::rowvec y(n-1,fill::zeros);
  for (int i=0;i<(n-1);i++){
    y(i) = x(i+1)-x(i);
  }
  return(y);
}
arma::rowvec Dtz(arma::rowvec z, const int N){
  arma::rowvec output(N,fill::zeros);
  arma::rowvec diffz = rcpp_diff(z);

  output(0) = -z(0);
  output(N-1) = z(N-2);
  for (int i=1;i<(N-1);i++){
    output(i) = -diffz(i-1);
  }
  return(output);
}

// 1. TVL2.IC
// [[Rcpp::export]]
arma::rowvec signal_tvl2_IC(arma::rowvec ytmp, const double lambda, const int maxiter){
  // 1-1. setup
  arma::rowvec y = ytmp;
  const int N = y.n_elem;
  arma::rowvec x(N,fill::zeros);
  arma::rowvec z(N-1,fill::zeros); // imaxiterialize
  const double alpha = 3.0;
  const double T = lambda/2;

  // 1-2. main iteration
  for (int k=0;k<maxiter;k++){
    x = y - Dtz(z, N);                // y-D'z
    z = z + (1.0/alpha)*rcpp_diff(x); // z + 1/alpha*D*z
    for (int i=0;i<(N-1);i++){
      if (z(i)>T){
        z(i) = T;
      }
      if (z(i)<-T){
        z(i) = -T;
      }
    }
  }
  // 1-3. return output
  return(x);
}

// 2. TVL2.MM
// [[Rcpp::export]]
arma::colvec signal_tvl2_MM(arma::colvec ytmp, const double lambda, const int maxiter){
  // 2-1. set up
  arma::colvec y = ytmp;
  const int N = y.n_elem;
  arma::mat D(N-1,N,fill::zeros);
  for (int i=0;i<(N-1);i++){
    D(i,i)   = -1;
    D(i,i+1) = 1;
  }
  arma::mat DDT = D*D.t();

  // 2-2. imaxiterialization
  arma::colvec x=y;
  arma::colvec Dx = D*x;
  arma::colvec Dy = D*y;
  arma::mat F(N-1,N-1,fill::zeros);

  // 2-3. main computation
  for (int k=0;k<maxiter;k++){
    F  = diagmat(abs(Dx)/lambda) + DDT;  // F : Sparse Banded Matrix
    x  = y - D.t()*solve(F,Dy);
    Dx = D*x;
  }

  // 2-4. return output
  return(x);
}
