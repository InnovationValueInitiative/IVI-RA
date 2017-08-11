// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// Convert armadillo vector to armadillo matrix by row
// [[Rcpp::export]]
arma::mat matrix_byrow(arma::rowvec v, int nrow, int ncol){
  int l = v.n_elem;
  arma::mat m1(0, l);
  m1.insert_rows(0, v);
  m1.reshape(ncol, nrow);
  m1 = arma::trans(m1);
  return(m1);
}

// Convert armadillo vector to armadillo matrix by column
// [[Rcpp::export]]
arma::mat matrix_bycol(arma::rowvec v, int nrow, int ncol){
  int l = v.n_elem;
  arma::mat m1(0, l);
  m1.insert_rows(0, v);
  m1.reshape(nrow, ncol);
  return(m1);
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
std::vector<std::vector<std::string> > charmat2stdvec(CharacterMatrix x){
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


// R list to arma::cube
// [[Rcpp::export]]
arma::cube array2cube(NumericVector array){
  IntegerVector dim = array.attr("dim"); 
  arma::cube my_array(array.begin(),dim[0], dim[1], dim[2], false);
  return my_array;
}