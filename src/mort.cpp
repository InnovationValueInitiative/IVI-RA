// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// Estimate probability given baseline probability and log of odds ratio
// [[Rcpp::export]]
double newprobC(arma::rowvec x, arma::rowvec logor, double logit_baseprob){
  double logit_factor = dot(x, logor);
  double newprob = 1/(1 + exp(-logit_baseprob - logit_factor));
  return newprob;
}

// Update mortality probability given change in baseline haq
// [[Rcpp::export]]
void update_qxC(double baseline_haq, double current_haq, double &qx, 
                       double cycle_length, double month, arma::rowvec loghr_vec){
  double change_haq = (current_haq - baseline_haq)/0.25;
  double rate = -log(1- qx);
  double loghr = 0;
  if (month <= 6){
    loghr = loghr_vec(0);
  }
  else if (month > 6 & month <=12){
    loghr = loghr_vec(1);
  }
  else if (month > 12 & month <= 24){
    loghr = loghr_vec(2);
  }
  else if (month > 24 & month <= 36){
    loghr = loghr_vec(3);
  }
  else if (month > 36){
    loghr = loghr_vec(4);
  }
  rate = rate * exp(loghr * change_haq);
  qx = 1 - exp(-rate * (cycle_length/12));
}

// Estimate mortality probability given patient characteristics at given point in time
// [[Rcpp::export]]
double mortprobC(int age, int male, arma::mat lifetable_male, arma::mat lifetable_female,
                 arma::rowvec x, arma::rowvec logor, double haq0, double haq,
                 double cycle_length, double month, arma::rowvec loghr_vec) {
  double logit_qxbase = 0;
  double qx = 0;
  
  // adjust probability of mortality using log odds ratio
  // then adjust for cycle length and change in baseline HAQ
  if (age >= 100){
    qx = 1;
  }
  else{
    // probability of mortality from lifetables
    if (male == 1){ 
      logit_qxbase = lifetable_male.row(age)(2);
    }
    else {
      logit_qxbase = lifetable_female.row(age)(2);
    }
    qx = newprobC(x, logor, logit_qxbase);
    update_qxC(haq0, haq, qx, cycle_length, month, loghr_vec);
  }
  
  return qx;
}

// Sample death (0/1 indicator) given patient characteristics at given point in time
// [[Rcpp::export]]
int sample_deathC(int age, int male, arma::mat lifetable_male, arma::mat lifetable_female,
                  arma::rowvec x, arma::rowvec logor, double haq0, double haq,
                  double cycle_length, double month, arma::rowvec loghr_vec) {
  double mortprob = mortprobC(age, male, lifetable_male, lifetable_female, x, logor,
                              haq0, haq, cycle_length, month, loghr_vec);
  unsigned int death = R::rbinom(1, mortprob);
  return death;
}

// Sample life-expectancy for a given age and a vector of HAQ scores using sample_deathC
// [[Rcpp::export]]
std::vector<double> sample_survC(int n, double age0, double male, arma::mat lifetable_male, arma::mat lifetable_female,
                  arma::rowvec x, arma::rowvec logor, double haq0, std::vector<double> haq,
                  double cycle_length, arma::rowvec loghr_vec) {
  int T = haq.size();
  std::vector<double> le;
  for (int i = 0; i < n; ++i){
    double age = age0;
    double month = 0;
    for (int t = 0; t < T; ++t){
      int d = sample_deathC(int(age), male, lifetable_male, lifetable_female,
                           x, logor, haq0, haq[t], cycle_length, month, loghr_vec);
      if (d == 1){
        break;
      }
      age = age + cycle_length/12;
      month = month + cycle_length;
    }
  le.push_back(age);
  }
 return le;
}

