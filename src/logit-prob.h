# ifndef LOGIT_PROB_H
# define LOGIT_PROB_H
#include <RcppArmadillo.h>

arma::rowvec ologit_probC(arma::rowvec x, arma::rowvec beta, arma::rowvec cut);
arma::rowvec mlogit_probC(arma::rowvec x, arma::mat beta);

# endif

