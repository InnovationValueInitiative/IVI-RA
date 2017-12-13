# ifndef LOGIT_PROB_H
# define LOGIT_PROB_H
#include <RcppArmadillo.h>

arma::rowvec mlogit_probC(arma::rowvec x, arma::mat beta);
double logistic(double p);

# endif

