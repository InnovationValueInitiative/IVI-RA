// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// Predicted probabilities for single individual from ordered logistic regression
//' @export
// [[Rcpp::export]]
arma::rowvec ologit_probC(arma::rowvec x, arma::rowvec beta, arma::rowvec cut) {
  int ncuts = cut.n_elem;
  double pstar = 0;
  double xb = dot(x, beta);
  double pstar_old = 0;
  arma::rowvec p(ncuts + 1);
      
  // loop over cut points
  for (int j = 0; j <= ncuts; ++j){
    if (j < ncuts){
      double invlogit = xb - arma::as_scalar(cut(j));
      pstar = 1/(1 + exp(invlogit));
      p(j) = pstar - pstar_old;
    }
    else{
      p(j) = 1 - pstar_old;
    }
    pstar_old = pstar;
  }
  return p;
}

// Predicted probabilities for single individual from multinomial logistic regression
//' @export
// [[Rcpp::export]]
arma::rowvec mlogit_probC(arma::rowvec x, arma::mat beta) {
  int ncat = beta.n_rows + 1;
  int zero = 0;
  arma::rowvec odds(ncat);
  arma::rowvec p;
  for (int j = 0; j < ncat; ++j){
    if (j == 0){
      odds(j) = exp(zero);
    }
    else{
      odds(j) = exp(dot(x, beta.row(j - 1)));
    }
  }
  double odds_sum = sum(odds);
  p = odds/odds_sum;
  return p;
}

// Logistic function
// [[Rcpp::export]]
double logistic(double p){
  return 1/(1 + exp(-p));
}

