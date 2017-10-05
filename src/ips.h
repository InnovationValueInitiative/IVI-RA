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

class ModelStructure {
public:
  std::string tx_ihaq;
  std::string tx_iswitch;
  std::string cdmards_haq_model;
  std::string ttd_cause;
  std::string ttd_dist;
  std::string utility_model;
  void set_model_structure(std::vector<std::string> x);
};

struct TTDPars1 {
  arma::rowvec loc;
  double anc1;
  double anc2;
};

class TTD{
private: 
  arma::rowvec x_eular;
  arma::rowvec x_all;
  arma::rowvec x_da;
  TTDPars1 eular_mod_pars;
  TTDPars1 eular_good_pars;
  TTDPars1 da_pars;
  TTDPars1 all_pars;
  int eular;
  int da_cat;
  int tswitch;
  double cycle_length;
  double ttsi;
  ModelStructure model_structure;
public:
  void set_x_ttd_da();
  double sim_ttd_eular();
  double sim_ttd_all();
  double sim_ttd_da();
  double sim_ttd();
  void set(arma::rowvec x_eular_, arma::rowvec x_all_, arma::rowvec x_da_,
           TTDPars1 eular_mod_pars_, TTDPars1 eular_good_pars_, 
           TTDPars1 da_pars_, TTDPars1 all_pars_, int eular_,
           int da_cat_, int tswitch_, double cycle_length_, 
           double ttsi_, ModelStructure model_structure_);
};

# endif
