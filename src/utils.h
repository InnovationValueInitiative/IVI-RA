# ifndef UTILS_H
# define UTILS_H
#include <Rcpp.h>

std::vector<std::vector<std::string> > stdvec2matrix(std::vector<std::string> x, int nrow, int ncol);
std::vector<std::vector<std::string> > charmat2stdvec(Rcpp::CharacterMatrix x);
double discount_factor(double t, double discount);
double annualized_value(double outcome, double df, int T);

# endif

