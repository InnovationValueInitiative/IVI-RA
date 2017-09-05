# ifndef ips_H
# define ips_H
#include <RcppArmadillo.h>

class nmaACR {
public:
  std::string hist;
  double rr;
  double A;
  double z2;
  double z3;
  arma::rowvec d_beta;
  arma::rowvec x;
  void set(std::string hist_, double rr_, double A_, double z2_, double z3_, 
           arma::rowvec d_beta_, arma::rowvec x_, int line);
  arma::rowvec nma_acrprob();
  double sim_acr();
};

class nmaLM {
public:
  std::string hist;
  double rr;
  double A;
  arma::rowvec d_beta;
  arma::rowvec x;
  void set(std::string hist_, double rr_, double A_, 
           arma::rowvec d_beta_, arma::rowvec x_, int line);
  double sim_dy();
};

class TxIHaq {
public:
  TxIHaq(): acr(NA_REAL), eular(NA_REAL), dhaq(0.0){ } 
  int acr;
  int eular;
  double dhaq;
  void sim(std::string tx_ihaq_model, int line, int therapy, int nbt,
           nmaACR nma_acr, nmaLM nma_dhaq,
           arma::mat acr2eular, arma::rowvec acr2haq, arma::rowvec eular2haq);
  
};

class TxISwitch {
public:
  int tswitch;
  double das28;
  double sdai;
  double cdai;
  int da_cat;
  int get_das28_cat();
  int get_sdai_cat();
  int get_cdai_cat();
  TxISwitch(): tswitch(0), das28(NA_REAL), sdai(NA_REAL), cdai(NA_REAL), da_cat(0){ } 
  double get_da_new(double da_old, double da_change, double lower, double upper);
  void sim(std::string model, int line, int therapy, int nbt, int acr, int eular,
           arma::rowvec acr2das28, arma::rowvec acr2sdai, arma::rowvec acr2cdai,
           nmaLM nma_das28, arma::rowvec p);
};


# endif
