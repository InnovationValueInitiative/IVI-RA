#include <Rcpp.h>

// Discount factor
double discount_factor(double t, double discount, double period_length = 1){
  double beta = 1/pow((1 + discount), period_length);
  return pow(beta, t);
}

// Calculate constant from geometric series
double const_geometric_series(double sum, double r, int T){
  return sum * (1 - r)/(1 - pow(r, T));
}

// Convert std::vector to matrix (ie., 2-dimensional vector) filling rowwise
std::vector<std::vector<std::string> > stdvec2matrix(std::vector<std::string> x, int nrow, int ncol){
  std::vector<std::vector<std::string> > vecs(nrow, std::vector<std::string>(ncol));
  for (int i = 0; i < ncol; ++i){
    for (int j = 0; j < nrow; ++j){
      vecs[i][j] = x[j + i * ncol];
    }
  }
  return vecs;
}

// Convert character matrix to 2-dimensional std string vector
std::vector<std::vector<std::string> > charmat2stdvec(Rcpp::CharacterMatrix x){
  int N = x.ncol();
  int M = x.nrow();
  std::vector<std::vector<std::string> > vecs(M, std::vector<std::string>(N));
  for (int i = 0; i < M; ++i){
    for (int j = 0; j < N; ++j){
      vecs[i][j] = x(i, j);
    }
  }
  return vecs;
}