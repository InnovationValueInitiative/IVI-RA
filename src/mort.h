# ifndef mort_H
# define mort_H
#include <RcppArmadillo.h>

double or2rrC(arma::rowvec x, arma::rowvec logor, double baseprob);

double newprobC(arma::rowvec x, arma::rowvec logor, double logit_baseprob);
  
double mortprobC(int age, int male, arma::mat lifetable_male, arma::mat lifetable_female,
                 arma::rowvec x, arma::rowvec logor);

void update_qxC(double baseline_haq, double current_haq, double qx, 
                double cycle_length, double month);

int sample_deathC(int age, int male, arma::mat lifetable_male, arma::mat lifetable_female,
                  arma::rowvec x, arma::rowvec logor, double haq0, double haq,
                  double cycle_length, double month, arma::rowvec loghr_vec);
# endif

