# ifndef RAND_H
# define RAND_H
#include <RcppArmadillo.h>

double rsurvC(double location, double anc1, std::string dist, double anc2 = 0.0);

double rbvcnormC(double &x, double &y_mean, double &x_mean, 
                 double &y_var, double &x_var, double &cor);

# endif

