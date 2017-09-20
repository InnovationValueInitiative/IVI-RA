#include <Rcpp.h>

// Discount factor
double discount_factor(double t, double discount){
  return 1/pow((1 + discount), t);
}

// Calculate annualized value
double annualized_value(double outcome, double df, int T){
  return outcome * (1 - df)/(1 - pow(df, T));
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