#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//[[Rcpp::export]]
NumericVector one_or_exp(NumericVector x) {
  NumericVector y(x);
  y = ifelse(x < 0.0, 1.0, exp(x));
  return y;
}

// [[Rcpp::export]]
List randomWalk2Rcpp(double niter, double lambda) {
  NumericVector x(niter);
  x[0] = 0;
  for (int i=1; i < niter; ++i) {
    x[i] = x[i-1] + lambda*(2.0 * Rf_rbinom(1, 0.5) - 1.0);
  }
  NumericVector y(niter);
  y[0] = 0;
  for (int i=1; i < niter; ++i) {
    y[i] = y[i-1] + lambda*(2.0 * Rf_rbinom(1, 0.5) - 1.0);
  }
  List z = List::create(Named("x") = x, _["y"] = y );
  return z ;
}

// [[Rcpp::export]]
arma::mat armadillo_solve(arma::mat A, arma::vec b) {
  int m = A.n_cols;
  arma::mat Z(m,1);
  Z = solve(A, b);
  return Z;
}

// [[Rcpp::export]]
arma::mat col_ridge_2(arma::mat Y, arma::mat X, arma::vec lambda) {
  int m = X.n_cols;
  int n = Y.n_cols;
  arma::mat ID = arma::eye<arma::mat>(m, m);
  arma::mat Betahat(m, n);
  for (int i=0; i < n; ++i) {
   Betahat.col(i) = (X.t() * X + lambda[i] * ID).i() * X.t() * Y.col(i);
  }
  return Betahat;
}

// [[Rcpp::export]]
arma::mat S_tilde_zero_Rcpp(arma::mat U, arma::vec trt, arma::vec time, arma::vec status, double n, arma::vec theta, double beta1, arma::vec beta2, double sigma, double t, arma::mat alpha, double maxm, double p) {
  arma::mat first = arma::zeros<arma::mat>(1,1);
  arma::mat Di(1,1);
  arma::mat Ui(10,p);
  NumericVector tnum=wrap(time>=t);
  for (int i=0; i < n; ++i) {
    arma::mat Ui = U.rows((i*maxm), ((i+1)*maxm-1));
    arma::mat Di = (theta + beta2*trt[i]).t() * (inv(Ui.t() * Ui)) * (theta + beta2*trt[i]);
    first += (tnum[i]) * (exp(beta1 * trt[i])) * exp(((theta+beta2*trt[i]).t() * ((alpha.row(i)).t())) - 0.5*sigma*sigma*Di);
  }
  return (first/n);
}


// [[Rcpp::export]]
arma::mat S_tilde_oneone_Rcpp(arma::mat U, arma::vec trt, arma::vec time, arma::vec status, double n, arma::vec theta, double beta1, arma::vec beta2, double sigma, double t, arma::mat alpha, double maxm, double p) {
  arma::mat first = arma::zeros<arma::mat>(p,1);
  arma::mat Di(1,1);
  arma::mat Ui(10,p);
  NumericVector tnum=wrap(time>=t);
  for (int i=0; i < n; ++i) {
    arma::mat Ui = U.rows((i*maxm), ((i+1)*maxm-1));
    arma::mat Di = (theta + beta2*trt[i]).t() * (inv(Ui.t() * Ui)) * (theta + beta2*trt[i]);
    first += as_scalar((tnum[i]) * (exp(beta1 * trt[i])) * exp(((theta+beta2*trt[i]).t() * ((alpha.row(i)).t())) - 0.5*sigma*sigma*Di)) * (((alpha.row(i)).t()) - sigma*sigma*(inv(Ui.t() * Ui))*(theta+beta2*trt[i]));
      }
  return (first/n);
}

// [[Rcpp::export]]
arma::mat S_tilde_onetwo_Rcpp(arma::mat U, arma::vec trt, arma::vec time, arma::vec status, double n, arma::vec theta, double beta1, arma::vec beta2, double sigma, double t, arma::mat alpha, double maxm, double p) {
  arma::mat first = arma::zeros<arma::mat>(1,1);
  arma::mat Di(1,1);
  arma::mat Ui(10,p);
  NumericVector tnum=wrap(time>=t);
  for (int i=0; i < n; ++i) {
    arma::mat Ui = U.rows((i*maxm), ((i+1)*maxm-1));
    arma::mat Di = (theta + beta2*trt[i]).t() * (inv(Ui.t() * Ui)) * (theta + beta2*trt[i]);
    first += as_scalar((tnum[i]) * (trt[i]) * (exp(beta1 * trt[i])) * exp(((theta+beta2*trt[i]).t() * ((alpha.row(i)).t())) - 0.5*sigma*sigma*Di));
  }
  return (first/n);
}

// [[Rcpp::export]]
arma::mat S_tilde_onethree_Rcpp(arma::mat U, arma::vec trt, arma::vec time, arma::vec status, double n, arma::vec theta, double beta1, arma::vec beta2, double sigma, double t, arma::mat alpha, double maxm, double p) {
  arma::mat first = arma::zeros<arma::mat>(p,1);
  arma::mat Di(1,1);
  arma::mat Ui(10,p);
  NumericVector tnum=wrap(time>=t);
  for (int i=0; i < n; ++i) {
    arma::mat Ui = U.rows((i*maxm), ((i+1)*maxm-1));
    arma::mat Di = (theta + beta2*trt[i]).t() * (inv(Ui.t() * Ui)) * (theta + beta2*trt[i]);
    first += as_scalar((tnum[i]) * (trt[i]) * (exp(beta1 * trt[i])) * exp(((theta+beta2*trt[i]).t() * ((alpha.row(i)).t())) - 0.5*sigma*sigma*Di)) * (((alpha.row(i)).t()) - sigma*sigma*(inv(Ui.t() * Ui))*(theta+beta2*trt[i]));
  }
  return (first/n);
}

// [[Rcpp::export]]
arma::mat S_tilde_zero_prime_Rcpp(arma::mat U, arma::vec trt, arma::vec time, arma::vec status, double n, arma::vec theta, double beta1, arma::vec beta2, double sigma, double t, arma::mat alpha, double maxm, double p) {
  arma::mat first = arma::zeros<arma::mat>((2*p+1),1);
  arma::mat convert = arma::ones<arma::mat>(1,1);
  arma::mat Di(1,1);
  arma::mat Ui(10,p);
  arma::mat comp1(p,1);
  arma::mat comp2(1,1);
  arma::mat comp3(p,1);
  NumericVector tnum=wrap(time>=t);
  for (int i=0; i < n; ++i) {
    arma::mat Ui = U.rows((i*maxm), ((i+1)*maxm-1));
    arma::mat Di = (theta + beta2*trt[i]).t() * (inv(Ui.t() * Ui)) * (theta + beta2*trt[i]);
    arma::mat comp1 = ((alpha.row(i)).t()) - sigma*sigma*(inv(Ui.t() * Ui))*(theta+beta2*trt[i]);
    arma::mat comp2 = (trt[i])*convert;
    arma::mat comp3 = (trt[i]) * ((alpha.row(i)).t()) - sigma*sigma*(inv(Ui.t() * Ui))*(theta+beta2*trt[i]);
    first += as_scalar((tnum[i]) * (exp(beta1 * trt[i])) * exp(((theta+beta2*trt[i]).t() * ((alpha.row(i)).t())) - 0.5*sigma*sigma*Di)) * join_vert(comp1, comp2, comp3);
  }
  return (first/n);
}

// [[Rcpp::export]]
arma::mat S_tilde_one_prime_one_Rcpp(arma::mat U, arma::vec trt, arma::vec time, arma::vec status, double n, arma::vec theta, double beta1, arma::vec beta2, double sigma, double t, arma::mat alpha, double maxm, double p) {
  arma::mat Di(1,1);
  arma::mat Ui(10,p);
  arma::mat comp1 = arma::zeros<arma::mat>(p,p);
  arma::mat comp2 = arma::zeros<arma::mat>(p,1);
  arma::mat comp3 = arma::zeros<arma::mat>(p,p);
  arma::mat first(p, (2*p+1));
  NumericVector tnum=wrap(time>=t);
  for (int i=0; i < n; ++i) {
    arma::mat Ui = U.rows((i*maxm), ((i+1)*maxm-1));
    arma::mat Di = (theta + beta2*trt[i]).t() * (inv(Ui.t() * Ui)) * (theta + beta2*trt[i]);
    comp1 += ((tnum[i])*(((alpha.row(i)).t()) - sigma*sigma*(inv(Ui.t() * Ui))*(theta+beta2*trt[i])) * ((((alpha.row(i)).t()) - sigma*sigma*(inv(Ui.t() * Ui))*(theta+beta2*trt[i])).t())* as_scalar((exp(beta1 * trt[i])) * exp(((theta+beta2*trt[i]).t() * ((alpha.row(i)).t())) - 0.5*sigma*sigma*Di)) - (tnum[i])*sigma*sigma*(inv(Ui.t() * Ui))*as_scalar((exp(beta1 * trt[i])) * exp(((theta+beta2*trt[i]).t() * ((alpha.row(i)).t())) - 0.5*sigma*sigma*Di)));
    comp2 += trt[i]*((tnum[i])*(((alpha.row(i)).t()) - sigma*sigma*(inv(Ui.t() * Ui))*(theta+beta2*trt[i])) * ((as_scalar((exp(beta1 * trt[i])) * exp(((theta+beta2*trt[i]).t() * ((alpha.row(i)).t())) - 0.5*sigma*sigma*Di)))));
    comp3 += trt[i]*((tnum[i])*(((alpha.row(i)).t()) - sigma*sigma*(inv(Ui.t() * Ui))*(theta+beta2*trt[i])) * ((((alpha.row(i)).t()) - sigma*sigma*(inv(Ui.t() * Ui))*(theta+beta2*trt[i])).t())* as_scalar((exp(beta1 * trt[i])) * exp(((theta+beta2*trt[i]).t() * ((alpha.row(i)).t())) - 0.5*sigma*sigma*Di)) - (tnum[i])*sigma*sigma*(inv(Ui.t() * Ui))*as_scalar((exp(beta1 * trt[i])) * exp(((theta+beta2*trt[i]).t() * ((alpha.row(i)).t())) - 0.5*sigma*sigma*Di)));
    first = join_horiz(comp1, comp2, comp3);
  }
  return (first/n);
}

// [[Rcpp::export]]
arma::mat S_tilde_one_prime_two_Rcpp(arma::mat U, arma::vec trt, arma::vec time, arma::vec status, double n, arma::vec theta, double beta1, arma::vec beta2, double sigma, double t, arma::mat alpha, double maxm, double p) {
  arma::mat first = arma::zeros<arma::mat>(1,(2*p+1));
  arma::mat convert = arma::ones<arma::mat>(1,1);
  arma::mat Di(1,1);
  arma::mat Ui(10,p);
  arma::mat comp1(p,1);
  arma::mat comp2(1,1);
  arma::mat comp3(p,1);
  NumericVector tnum=wrap(time>=t);
  for (int i=0; i < n; ++i) {
    arma::mat Ui = U.rows((i*maxm), ((i+1)*maxm-1));
    arma::mat Di = (theta + beta2*trt[i]).t() * (inv(Ui.t() * Ui)) * (theta + beta2*trt[i]);
    arma::mat comp1 = ((alpha.row(i)).t()) - sigma*sigma*(inv(Ui.t() * Ui))*(theta+beta2*trt[i]);
    arma::mat comp2 = (trt[i])*convert;
    arma::mat comp3 = (trt[i]) * ((alpha.row(i)).t()) - sigma*sigma*(inv(Ui.t() * Ui))*(theta+beta2*trt[i]);
    first += as_scalar(trt[i]*(tnum[i]) * (exp(beta1 * trt[i])) * exp(((theta+beta2*trt[i]).t() * ((alpha.row(i)).t())) - 0.5*sigma*sigma*Di)) * (join_vert(comp1, comp2, comp3)).t();
  }
  return (first/n);
}

// [[Rcpp::export]]
arma::mat S_tilde_one_prime_three_Rcpp(arma::mat U, arma::vec trt, arma::vec time, arma::vec status, double n, arma::vec theta, double beta1, arma::vec beta2, double sigma, double t, arma::mat alpha, double maxm, double p) {
  arma::mat Di(1,1);
  arma::mat Ui(10,p);
  arma::mat comp1 = arma::zeros<arma::mat>(p,p);
  arma::mat comp2 = arma::zeros<arma::mat>(p,1);
  arma::mat comp3 = arma::zeros<arma::mat>(p,p);
  arma::mat first(p, (2*p+1));
  NumericVector tnum=wrap(time>=t);
  for (int i=0; i < n; ++i) {
    arma::mat Ui = U.rows((i*maxm), ((i+1)*maxm-1));
    arma::mat Di = (theta + beta2*trt[i]).t() * (inv(Ui.t() * Ui)) * (theta + beta2*trt[i]);
    comp1 += trt[i]*((tnum[i])*(((alpha.row(i)).t()) - sigma*sigma*(inv(Ui.t() * Ui))*(theta+beta2*trt[i])) * ((((alpha.row(i)).t()) - sigma*sigma*(inv(Ui.t() * Ui))*(theta+beta2*trt[i])).t())* as_scalar((exp(beta1 * trt[i])) * exp(((theta+beta2*trt[i]).t() * ((alpha.row(i)).t())) - 0.5*sigma*sigma*Di)) - (tnum[i])*sigma*sigma*(inv(Ui.t() * Ui))*as_scalar((exp(beta1 * trt[i])) * exp(((theta+beta2*trt[i]).t() * ((alpha.row(i)).t())) - 0.5*sigma*sigma*Di)));
    comp2 += trt[i]*((tnum[i])*(((alpha.row(i)).t()) - sigma*sigma*(inv(Ui.t() * Ui))*(theta+beta2*trt[i])) * ((as_scalar((exp(beta1 * trt[i])) * exp(((theta+beta2*trt[i]).t() * ((alpha.row(i)).t())) - 0.5*sigma*sigma*Di)))));
    comp3 += trt[i]*((tnum[i])*(((alpha.row(i)).t()) - sigma*sigma*(inv(Ui.t() * Ui))*(theta+beta2*trt[i])) * ((((alpha.row(i)).t()) - sigma*sigma*(inv(Ui.t() * Ui))*(theta+beta2*trt[i])).t())* as_scalar((exp(beta1 * trt[i])) * exp(((theta+beta2*trt[i]).t() * ((alpha.row(i)).t())) - 0.5*sigma*sigma*Di)) - (tnum[i])*sigma*sigma*(inv(Ui.t() * Ui))*as_scalar((exp(beta1 * trt[i])) * exp(((theta+beta2*trt[i]).t() * ((alpha.row(i)).t())) - 0.5*sigma*sigma*Di)));
    first = join_horiz(comp1, comp2, comp3);
  }
  return (first/n);
}







