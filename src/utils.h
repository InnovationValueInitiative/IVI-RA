# ifndef UTILS_H
# define UTILS_H
#include <RcppArmadillo.h>


std::vector<std::vector<std::string> > stdvec2matrix(std::vector<std::string> x, int nrow, int ncol);
std::vector<std::vector<std::string> > charmat2stdvec(Rcpp::CharacterMatrix x);
double discount_factor(double discount, double t);

# endif

