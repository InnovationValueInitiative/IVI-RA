#include <RcppArmadillo.h>
#include "ips.h"
using namespace Rcpp;

/*********************
* Testing nmaACR class
**********************/
// [[Rcpp::export]]
int sim_acr_test(){
  nmaACR sim;
  double A = 0;
  double z2 = .64;
  double z3 = 1.22;
  arma::rowvec d_beta(2); d_beta(0) = -.05; d_beta(1) = .05;
  arma::rowvec x(2); x(0) = 1; x(1) = 1;
  
  // (1) hist = naive
  std::string hist = "naive";
  double rr = .8;
  
  // line =  0
  int line = 0;
  sim.set(hist, rr, A, z2, z3, d_beta, x, line);
  arma::rowvec acr_prob = sim.acrprobs();
  if (acr_prob(0) != 0.5) {
    Rcpp::stop("Fail hist = naive, rr = .8, line = 0");
  }
  
  // line = 1
  line = 1;
  sim.set(hist, rr, A, z2, z3, d_beta, x, line);
  acr_prob = sim.acrprobs();
  if (acr_prob(0) != 0.6) {
    Rcpp::stop("Fail hist = naive, rr = .8, line = 1");
  }
  
  // rr = 0
  rr = 0;
  sim.set(hist, rr, A, z2, z3, d_beta, x, line);
  acr_prob = sim.acrprobs();
  if (acr_prob(0) != 1.0) {
    Rcpp::stop("Fail hist = naive, rr = 0, line = 1");
  }
  
  // (2) hist = experienced
  hist = "experienced";
  rr = .8;
  line = 0;
  sim.set(hist, rr, A, z2, z3, d_beta, x, line);
  acr_prob = sim.acrprobs();
  if (acr_prob(0) != 0.6) {
    Rcpp::stop("Fail hist = experienced, rr = .8, line = 0");
  }
  
  return 0;
}

/*********************
* Testing nmaLM class
**********************/
// [[Rcpp::export]]
int sim_lm_test(){
  nmaLM sim;
  double A = -.25;
  arma::rowvec d_beta(2); d_beta(0) = -.05; d_beta(1) = .05;
  arma::rowvec x(2); x(0) = 1; x(1) = 1;
  
  // (1) hist = naive
  std::string hist = "naive";
  double rr = .8;
  
  // line =  0
  int line = 0;
  sim.set(hist, rr, A, d_beta, x, line);
  if (sim.sim_dy() != -.25) {
    Rcpp::stop("Fail hist = naive, rr = .8, line = 0");
  }
  
  // line = 1
  line = 1;
  sim.set(hist, rr, A, d_beta, x, line);
  if (sim.sim_dy() != -.2) {
    Rcpp::stop("Fail hist = naive, rr = .8, line = 1");
  }
  
  // rr = 0
  rr = 0;
  sim.set(hist, rr, A, d_beta, x, line);
  if (sim.sim_dy() != 0) {
    Rcpp::stop("Fail hist = naive, rr = 0, line = 1");
  }
  
  // (2) hist = experienced
  hist = "experienced";
  rr = .8;
  line = 0;
  sim.set(hist, rr, A, d_beta, x, line);
  if (sim.sim_dy() != -.2) {
    Rcpp::stop("Fail hist = experienced, rr = .8, line = 0");
  }

  return 0;
}