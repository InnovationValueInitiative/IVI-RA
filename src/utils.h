# ifndef UTILS_H
# define UTILS_H
#include <RcppArmadillo.h>

std::vector<std::vector<std::string> > stdvec2matrix(std::vector<std::string> x, int nrow, int ncol);
std::vector<std::vector<std::string> > charmat2stdvec(Rcpp::CharacterMatrix x);
double discount_factor(double t, double discount, double period_length = 1);
double const_geometric_series(double sum, double r, int T, int start = 1);
bool arma_rowvec_anyNA(arma::rowvec x);

# endif

