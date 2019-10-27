#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

/*
 * Image 1 : TV-L2 using Finite Difference
 */
// [[Rcpp::export]]
arma::mat image_tvl2_FD(arma::mat u0tmp, const double lambda, const double niter){
  // 1. initialize and setup
  arma::mat u0 = u0tmp;
  const int M = u0.n_rows;
  const int N = u0.n_cols;
  const double h = 1.0;         // space discretization
  const double eps = 0.00001;   // regularize TV at the origin
  const double eps2= eps*eps;
  arma::mat u(M,N,fill::zeros); // 'u' as new denoised image

  // 2. Begin Iterations in ITER
  double ux, uy, Gradu, co, co1, co2, co3, co4, div;
  for (int iter=0;iter<niter;iter++){
    // 2-1. looping over interior domain
    for (int i=1;i<(M-1);i++){
      for (int j=1;j<(N-1);j++){
        // 2-1-1. computing coefficients co1, co2, co3, co4
        ux=(u(i+1,j)-u(i,j))/h;
        uy=(u(i,j+1)-u(i,j))/h;
        Gradu=sqrt(eps2+ux*ux+uy*uy);
        co1=1.0/Gradu;

        ux=(u(i,j)-u(i-1,j))/h;
        uy=(u(i-1,j+1)-u(i-1,j))/h;
        Gradu=sqrt(eps2+ux*ux+uy*uy);
        co2=1.0/Gradu;

        ux=(u(i+1,j)-u(i,j))/h;
        uy=(u(i,j+1)-u(i,j))/h;
        Gradu=sqrt(eps2+ux*ux+uy*uy);
        co3=1.0/Gradu;

        ux=(u(i+1,j-1)-u(i,j-1))/h;
        uy=(u(i,j)-u(i,j-1))/h;
        Gradu=sqrt(eps2+ux*ux+uy*uy);
        co4=1.0/Gradu;

        co=1.0+(1/(2*lambda*h*h))*(co1+co2+co3+co4);
        div=co1*u(i+1,j)+co2*u(i-1,j)+co3*u(i,j+1)+co4*u(i,j-1);
        u(i,j)=(1.0/co)*(u0(i,j)+(1/(2*lambda*h*h))*div);
      }
    }
    // 2-2. free boundary conditions in u
    for (int i=1;i<(M-1);i++){
      u(i,0)  =u(i,1);
      u(i,N-1)=u(i,N-2);
    }
    for (int j=1;j<(N-1);j++){
      u(0,j)  = u(1,j);
      u(M-1,j)= u(M-2,j);
    }
    u(0,0) = u(1,1);
    u(0,N-1) = u(1,N-2);
    u(M-1,0) = u(M-2,1);
    u(M-1,N-1) = u(M-2,N-2);
  }

  // return output
  return(u);
}


arma::mat rcpp_imgradient_hor(arma::mat u){
  const int height = u.n_rows;
  const int width  = u.n_cols;
  arma::mat gradu(height,width,fill::zeros);

  for (int i=0;i<(width-1);i++){
    gradu.col(i) = u.col(i+1);
  }
  gradu.col(width-1) = u.col(width-1);
  gradu = gradu - u;
  return(gradu);
}
arma::mat rcpp_imgradient_ver(arma::mat u){
  const int height = u.n_rows;
  const int width  = u.n_cols;
  arma::mat gradu(height,width,fill::zeros);

  for (int i=0;i<(height-1);i++){
    gradu.row(i) = u.row(i+1);
  }
  gradu.row(height-1) = u.row(height-1);
  gradu = gradu - u;
  return(gradu);
}
arma::mat rcpp_pmax_matrix(arma::mat u, double x){
  const int m = u.n_rows;
  const int n = u.n_cols;
  arma::mat output(m,n,fill::zeros);
  double tgt = 0.0;

  for (int i=0;i<m;i++){
    for (int j=0;j<n;j++){
      tgt = u(i,j);
      if (x>=tgt){
        output(i,j) = x;
      } else {
        output(i,j) = tgt;
      }
    }
  }
  return(output);
}

/*
* Image 2. TV-L1 Primal Dual
*/
// [[Rcpp::export]]
arma::mat image_tvl1_primaldual(arma::mat utmp, const double lambda, const double niter){
  // 1. set parameters
  arma::mat u = utmp;
  double L2 = 8.0;
  double tau = 0.02;
  double sigma = 1.0/(L2*tau);
  double theta = 1.0;
  double lt = lambda*tau;

  const int height = u.n_rows;
  const int width  = u.n_cols;

  // 2. set empty arrays
  arma::mat nim = u;
  arma::mat unew(height,width,fill::zeros);
  arma::mat v(height,width,fill::zeros);
  arma::mat p1(height,width,fill::zeros);
  arma::mat p2(height,width,fill::zeros);
  arma::mat d(height,width,fill::zeros);
  arma::mat ux(height,width,fill::zeros);
  arma::mat uy(height,width,fill::zeros);

  // 3. approximate values
  p1 = rcpp_imgradient_hor(u);
  p2 = rcpp_imgradient_ver(u);

  // 4. main iteration
  arma::mat tester(height, width);
  arma::mat normep(height, width);
  arma::mat div;
  arma::mat zeros1width(1, width, fill::zeros);
  arma::mat zerosheight1(height, 1, fill::zeros);
  for (int k=0;k<niter;k++){
    // 4-1. compute gradient
    ux = rcpp_imgradient_hor(u);
    uy = rcpp_imgradient_ver(u);
    p1 = p1 + sigma*ux;
    p2 = p2 + sigma*uy;

    // 4-2. project
    tester = sqrt(p1%p1 + p2%p2);
    normep = rcpp_pmax_matrix(tester, 1.0);
    p1 = p1/normep;
    p2 = p2/normep;

    // 4-3. shrinkage
    div = join_cols(p2.head_rows(height-1),zeros1width) - join_cols(zeros1width, p2.head_rows(height-1));
    div = join_rows(p1.head_cols(width-1),zerosheight1) - join_rows(zerosheight1, p1.head_cols(width-1)) + div;

    // 4-4. TV-L1 model
    v = u + tau*div;
    for (int i=0;i<height;i++){
      for (int j=0;j<width;j++){
        // 3 term, seperately
        if (v(i,j)-nim(i,j) > lt){
          unew(i,j) = v(i,j)-lt;
        } else {
          unew(i,j) = 0;
        }
        if (v(i,j)-nim(i,j) < -lt){
          unew(i,j) += v(i,j)+lt;
        }
        if (pow(v(i,j)-nim(i,j),2)<=pow(lt,2)){
          unew(i,j) += nim(i,j);
        }
      }
    }

    // 4-5. extragradient step
    u = unew + theta*(unew-u);
  }
  // 5. return output matrix
  return(u);
}


/*
* Image 3. TV-L2 Primal Dual
*/
// [[Rcpp::export]]
arma::mat image_tvl2_primaldual(arma::mat utmp, const double lambda, const double niter){
  // 1. set parameters
  arma::mat u = utmp;
  double L2 = 8.0;
  double tau = 0.02;
  double sigma = 1.0/(L2*tau);
  double theta = 1.0;
  double lt = lambda*tau;

  const int height = u.n_rows;
  const int width  = u.n_cols;

  // 2. set empty arrays
  arma::mat nim = u;
  arma::mat unew(height,width,fill::zeros);
  arma::mat v(height,width,fill::zeros);
  arma::mat p1(height,width,fill::zeros);
  arma::mat p2(height,width,fill::zeros);
  arma::mat d(height,width,fill::zeros);
  arma::mat ux(height,width,fill::zeros);
  arma::mat uy(height,width,fill::zeros);

  // 3. approximate values
  p1 = rcpp_imgradient_hor(u);
  p2 = rcpp_imgradient_ver(u);

  // 4. main iteration
  arma::mat tester(height, width);
  arma::mat normep(height, width);
  arma::mat div;
  arma::mat zeros1width(1, width, fill::zeros);
  arma::mat zerosheight1(height, 1, fill::zeros);
  for (int k=0;k<niter;k++){
    // 4-1. compute gradient
    ux = rcpp_imgradient_hor(u);
    uy = rcpp_imgradient_ver(u);
    p1 = p1 + sigma*ux;
    p2 = p2 + sigma*uy;

    // 4-2. project
    tester = sqrt(p1%p1 + p2%p2);
    normep = rcpp_pmax_matrix(tester, 1.0);
    p1 = p1/normep;
    p2 = p2/normep;

    // 4-3. shrinkage
    div = join_cols(p2.head_rows(height-1),zeros1width) - join_cols(zeros1width, p2.head_rows(height-1));
    div = join_rows(p1.head_cols(width-1),zerosheight1) - join_rows(zerosheight1, p1.head_cols(width-1)) + div;

    // 4-4. TV-L2 model
    unew = (u + tau*div + lt*nim)/(1+tau);

    // 4-5. extragradient step
    u = unew + theta*(unew-u);
  }
  // 5. return output matrix
  return(u);
}
