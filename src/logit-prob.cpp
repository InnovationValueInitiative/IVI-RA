// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// Predicted probabilities for single individual from multinomial logistic regression
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


