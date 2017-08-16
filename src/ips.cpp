// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(hesim)]]

#include <RcppArmadillo.h>
#include "rand.h"
#include "mort.h"
#include "logit-prob.h"
#include "utils.h"
#include <hesim.h>
using namespace Rcpp;

// Model structure
struct ModelStructure {
  std::string tx_ihaq;
  std::string tx_iswitch;
  std::string cdmards_haq_model;
  std::string ttd_dist;
  std::string utility_model;
};

ModelStructure update_model_structure(ModelStructure m, std::vector<std::string> x){
  m.tx_ihaq = x[0];
  m.tx_iswitch = x[1];
  m.cdmards_haq_model = x[2];
  m.ttd_dist = x[3];
  m.utility_model = x[4];
  return m;
}

// Update HAQ initial response
void update_haq_t1(double &haq, double haq_change){
  haq = haq + haq_change;
}

// Effect of treatment on HAQ during the initial treatment phase
struct TxIHaq {
  TxIHaq(): acr(NA_REAL), eular(NA_REAL), dhaq(0.0){ } 
  int acr;
  int eular;
  double dhaq;
};

TxIHaq sim_tx_ihaq(std::string tx_ihaq_model, int line, int therapy, int nbt,
                     arma::rowvec nma_acr1, arma::rowvec nma_acr2, 
                     double nma_dhaq1, double nma_dhaq2,
                     arma::mat acr2eular, arma::rowvec acr2haq, arma::rowvec eular2haq){
  TxIHaq sim;
  
  // ACR response
  if (tx_ihaq_model == "acr-haq" || tx_ihaq_model == "acr-eular-haq"){
      if (therapy != nbt){
        if (line == 0){
          sim.acr = hesim::rcat1C(nma_acr1);
        }  
        else{
          sim.acr = hesim::rcat1C(nma_acr2);
        }
      }
      else{
          sim.acr = 0;
      }
  }
  
  // EULAR response
  if (tx_ihaq_model == "acr-eular-haq"){
      if (therapy != nbt){
        sim.eular = hesim::rcat1C(acr2eular.row(sim.acr));
      }
      else{
        sim.eular = 0;
      }
  }
  
  // HAQ
  if (tx_ihaq_model == "acr-haq" && therapy != nbt){
        sim.dhaq = acr2haq(sim.acr);
  } 
  else if (tx_ihaq_model == "acr-eular-haq" && therapy != nbt){
      sim.dhaq = eular2haq(sim.eular);
  } 
  else if (tx_ihaq_model == "haq" && therapy != nbt) {
      if (line == 0){
        sim.dhaq = nma_dhaq1;
      } else {
        sim.dhaq = nma_dhaq2;
      }
  }
  else if (therapy == nbt){
      sim.dhaq = 0.0;
  }
  
  // return
  return sim;
}

// [[Rcpp::export]]
List test_tx_ihaq(std::string tx_ihaq_model, int line, int therapy, int nbt,
                      arma::rowvec nma_acr1, arma::rowvec nma_acr2, 
                      double nma_dhaq1, double nma_dhaq2,
                      arma::mat acr2eular, arma::rowvec acr2haq, arma::rowvec eular2haq){
  TxIHaq sim = sim_tx_ihaq(tx_ihaq_model, line, therapy, nbt,
                        nma_acr1, nma_acr2, nma_dhaq1, nma_dhaq2,
                        acr2eular, acr2haq, eular2haq);
  return List::create(Rcpp::Named("acr") = sim.acr, Rcpp::Named("eular") = sim.eular,
                      Rcpp::Named("dhaq") = sim.dhaq);
}

// Update HAQ after initial response with linear progression
void update_haq_t(double &haq, double haq_change_therapy, 
                  arma::rowvec haq_change_age_vec, double age, double cycle_length){
  double haq_change_age = 0.0;
  if (age < 40){
    haq_change_age = haq_change_age_vec(0);
  }
  else if (age >= 40 & age <= 64){
    haq_change_age = haq_change_age_vec(1);
  }
  else {
    haq_change_age = haq_change_age_vec(2);
  }
  haq = haq + (haq_change_age + haq_change_therapy) * (cycle_length/12);
}

// Effect of treatment on switching treatment during the treatment phase
struct Tx_ISwitch {
  Tx_ISwitch(): tswitch(0), das28(NA_REAL), sdai(NA_REAL), cdai(NA_REAL), da_cat(0){ } 
  int tswitch;
  double das28;
  double sdai;
  double cdai;
  int da_cat;
};

// [[Rcpp::export]]
int get_das28_cat(double das28){
  int cat = 0;
  if (das28 >= 2.6 & das28 < 3.2){
    cat = 1;
  }
  else if (das28 >= 3.2 & das28 <= 5.1){
    cat = 2;
  }
  else if (das28 > 5.1){
    cat = 3;
  }
  return cat;
}

// [[Rcpp::export]]
int get_sdai_cat(double sdai){
  int cat = 0;
  if (sdai > 3.3 & sdai <= 11.0){
    cat = 1;
  }
  else if (sdai > 11 & sdai <= 26){
    cat = 2;
  }
  else if (sdai > 26){
    cat = 3;
  }
  return cat;
}

// [[Rcpp::export]]
int get_cdai_cat(double cdai){
  int cat = 0;
  if (cdai > 2.8 & cdai <= 10.0){
    cat = 1;
  }
  else if (cdai > 10 & cdai <= 22){
    cat = 2;
  }
  else if (cdai > 22){
    cat = 3;
  }
  return cat;
}

// [[Rcpp::export]]
double get_da_new(double da_old, double da_change, double lower, double upper){
  double da_new = da_old + da_change;
  if (da_new < lower){
    da_new = 0;
  } else if (da_new >upper){
    da_new = upper;
  }
  return da_new;
}

Tx_ISwitch sim_tx_iswitch(std::string tx_iswitch_model, int line, int therapy, int nbt,
                  int acr, int eular, double das28, double sdai, double cdai,
                  arma::rowvec acr2das28, arma::rowvec acr2sdai,arma::rowvec acr2cdai,
                  arma::rowvec nma_das28_1, arma::rowvec nma_das28_2, arma::rowvec p){
  
  Tx_ISwitch sim;
  sim.tswitch = 0;
  sim.da_cat = 0;
  double da_change = 0.0;
  
  // Switching rules
  if (tx_iswitch_model == "acr-switch"){ // Treatment -> ACR -> Switch
    if (acr == 0){
      sim.tswitch = 1;
    }
  }
  else if (tx_iswitch_model == "acr-das28-switch" ||   // Treatment -> ACR -> DA -> Switch
           tx_iswitch_model == "acr-sdai-switch" ||
           tx_iswitch_model == "acr-cdai-switch") {
      if (tx_iswitch_model == "acr-das28-switch"){
        da_change = acr2das28(acr);
        sim.das28 = get_da_new(das28, da_change, 0, 9.4);
        sim.da_cat = get_das28_cat(sim.das28);
      }
      else if (tx_iswitch_model == "acr-sdai-switch"){
        da_change = acr2sdai(acr);
        sim.sdai = get_da_new(sdai, da_change, 0, 86);
        sim.da_cat = get_sdai_cat(sdai);
      }
      else if (tx_iswitch_model == "acr-cdai-switch"){
        da_change = acr2cdai(acr);
        sim.cdai = get_da_new(cdai, da_change, 0, 76);
        sim.da_cat = get_cdai_cat(sim.cdai);
      }
      sim.tswitch = R::rbinom(1, p(sim.da_cat));
    }
  else if (tx_iswitch_model == "das28-switch"){ // Treatment -> DA -> Switch
      if (line == 0){
        sim.das28 = das28 + nma_das28_1(therapy);
        sim.da_cat = get_das28_cat(sim.das28);
      }
      else{
        sim.das28 = das28 + nma_das28_2(therapy);
        sim.da_cat = get_das28_cat(sim.das28);
      }
      sim.tswitch = R::rbinom(1, p(sim.da_cat));
  }
  else if (tx_iswitch_model == "acr-eular-switch"){ // Treatment -> ACR -> EULAR -> Switch
    if (eular == 0){
      sim.tswitch = 1;
    }
  }
  
  // NBT
  if (therapy == nbt){
    sim.tswitch = 0;
  }
  return sim;
}

//' @export
// [[Rcpp::export]]
List test_tx_iswitch(std::string tx_iswitch_model, int line, int therapy, int nbt,
                        int acr, int eular, double das28, double sdai, double cdai,
                        arma::rowvec acr2das28, arma::rowvec acr2sdai, arma::rowvec acr2cdai,
                        arma::rowvec nma_das28_1, arma::rowvec nma_das28_2, arma::rowvec p){
  Tx_ISwitch sim = sim_tx_iswitch(tx_iswitch_model, line, therapy, nbt, acr, eular,
                                     das28, sdai, cdai, acr2das28, acr2sdai, acr2cdai,
                                     nma_das28_1, nma_das28_2, p);
  return List::create(Rcpp::Named("tswitch") = sim.tswitch, Rcpp::Named("das28") = sim.das28,
                      Rcpp::Named("sdai") = sim.sdai, Rcpp::Named("cdai") = sim.cdai);
}


// Time to treatment discontinuation by eular response 
//' @export
// [[Rcpp::export]]
double sim_ttd_eular(arma::rowvec x, arma::rowvec loc_mod, double anc1_mod,
                 arma::rowvec loc_good, double anc1_good, 
                 int eular, std::string dist, double cycle_length, double ttsi,
                 double anc2_mod = 0.0, double anc2_good = 0.0){
  double surv = 0.0;
  if (ttsi < 0){
    surv = 0.0;
  }
  else {
    if (eular == 0){ //no eular response
      surv = 0.0;
    }
    else if (eular == 1){ //moderate eular responder
      surv = rsurvC(dot(x, loc_mod), anc1_mod, dist, anc2_mod);
    }
    else if (eular == 2){ //good eular responder
      surv = rsurvC(dot(x, loc_good), anc1_good, dist, anc2_good);
    }
  }
  return surv/cycle_length; // surv is measured in years, so surv/cycle_length is measured in model cycles
}

// Time to treatment discontinuation single model
//' @export
// [[Rcpp::export]]
double sim_ttd(arma::rowvec x, arma::rowvec loc, double anc1,
                   int tswitch, std::string dist,
                  double cycle_length, double ttsi,
                  double anc2 = 0.0, int da_cat = 0){
  double surv = 0.0;
  if (ttsi < 0 || tswitch == 1){
    surv = 0.0;
  }
  else {
      surv = rsurvC(dot(x, loc), anc1, dist, anc2);
  }
  return surv/cycle_length; // surv is measured in years, so surv/cycle_length is measured in model cycles
}

arma::rowvec update_x_ttd_da(arma::rowvec x, int da_cat){
  if (da_cat == 3){
    x(2) = 1;
  }
  else{
    x(2) = 0;
  }
  if (da_cat == 2){
    x(1) = 1;
  }
  else{
    x(1) = 0;
  }
  return x;
}

// Convert R TTD parameters to vectors of armadilo objects
struct TTDPars {
  std::map<std::string, arma::mat> loc;
  std::map<std::string, arma::vec> anc1;
  std::map<std::string, arma::vec> anc2;
};

TTDPars get_ttd_pars(Rcpp::List x){
  TTDPars ttd;
  int N = x.size();
  std::vector<std::string> dist = x.attr("names");
  for (int i = 0; i < N; ++i){
    Rcpp::List x_dist = x[i];
    arma::mat sample = as<arma::mat>(x_dist["sample"]);
    int N = sample.n_rows;
    arma::uvec loc_indx = as<arma::uvec> (x_dist["loc.index"]);
    int anc1_indx = as<int> (x_dist["anc1.index"]);
    int anc2_indx = as<int> (x_dist["anc2.index"]);
    arma::mat loc = sample.cols(loc_indx - 1);
    arma::vec anc1(N);
    arma::vec anc2(N);
    if(anc1_indx != NA_INTEGER){
      anc1 = sample.col(anc1_indx -1);
    }
    if(anc2_indx != NA_INTEGER){
      anc2 = sample.col(anc2_indx - 1);
    }
    ttd.loc.insert(std::make_pair(dist[i], loc));
    ttd.anc1.insert(std::make_pair(dist[i], anc1));
    ttd.anc2.insert(std::make_pair(dist[i], anc2));
  }
  return ttd;
}

// [[Rcpp::export]]
Rcpp::List test_ttd(Rcpp::List x){
  TTDPars ttd = get_ttd_pars(x);
  return Rcpp::List::create(Rcpp::Named("loc") = ttd.loc,
                            Rcpp::Named("anc1") = ttd.anc1,
                            Rcpp::Named("anc2") = ttd.anc2);
}


// Simulate latent class from multinomial logistic regression
//' @export
// [[Rcpp::export]]
int sim_mlogit_classC(arma::rowvec w, arma::mat delta){
  arma::rowvec latclass_prob = mlogit_probC(w, delta);
  int latclass = hesim::rcat1C(latclass_prob); 
  return(latclass);
}

// Simulate change in HAQ from Norton mixture model conditional on class
//' @export
// [[Rcpp::export]]
double sim_dhaq_class_lcgm1C(double year, double cycle_length, arma::rowvec beta){
  double xt = 1 - (1/(1 + year));
  double year_lag = year - cycle_length/12;
  double xt_lag = 1 - (1/(1 + year_lag));
  double dyhat = beta(1) * (xt - xt_lag) + beta(2) * (pow(xt, 2) - pow(xt_lag, 2)) +
    beta(3) * (pow(xt, 3) - pow(xt_lag, 3));
  return(dyhat);
}

// Simulate change in HAQ from Norton mixture model 
//' @export
// [[Rcpp::export]]
double sim_dhaq_lcgm1C(double year, double cycle_length, double age, double female,
                       double das28, arma::mat delta, arma::mat beta){
  arma::rowvec w(8);
  w(0) = 1.0; w(1) = age; w(2) = female; w(3) = das28;
  w(4) = 8.2; // disease duration
  w(5) = .73; // rheumatoid factor
  w(6) = 1; // ACR 1987 criteria
  w(7) = .49; // IMDQ4 measure of SES
  int latclass = sim_mlogit_classC(w, delta);
  arma::rowvec beta_vec = beta.row(latclass);
  double dyhat = sim_dhaq_class_lcgm1C(year, cycle_length, beta_vec);
  return(dyhat);
}

// Simulate costs
struct Cost {
  double tx;
  double hosp;
  double mgmt;
  double si;
  double prod;
};

//// Treatment Costs
// [[Rcpp::export]]
double sim_tx_cost1C(int t, arma::rowvec agents_ind, std::vector<std::string> tx_name,
                     std::vector<double> init_dose_val, std::vector<double> ann_dose_val,
                     std::vector<double> strength_val,
                     std::vector<double> init_num_doses, std::vector<double> ann_num_doses,
                     std::vector<double> price, 
                     std::vector<double> infusion_cost, std::vector<int> loading_dose,
                     std::vector<int> weight_based, double weight, double cycle_length,
                     arma::rowvec discount){
  
  // initial vs. maintenance phase dosing
  std::vector<double> dose_amt;
  std::vector<double> dose_num;
  if (t == 0){
    dose_amt = init_dose_val;
    dose_num = init_num_doses;
  }
  else{
    dose_amt = ann_dose_val;
    dose_num = ann_num_doses;
  }

  // Calculate costs
  int J = agents_ind.size();
  double tc = 0;
  double infusion_cost_j = 0;
  for (int j = 0; j < J; ++j){
    if(arma::is_finite(agents_ind(j))){
      int agents_ind_j = agents_ind(j);
      double dose_num_j = dose_num[agents_ind_j];
      std::string tx_name_j = tx_name[agents_ind_j];
      if (tx_name_j == "tcz" && weight >= 100){
        dose_num_j = dose_num_j * 2;
      }
      if (weight_based[agents_ind_j] == 0){
        weight = 1; 
      }
      if (t == 0){
        if (loading_dose[agents_ind_j] == 1){
          infusion_cost_j = infusion_cost[agents_ind_j];
        } 
        else{
          infusion_cost_j = infusion_cost[agents_ind_j] * dose_num[agents_ind_j];
        }
        tc +=  ceil(weight * dose_amt[agents_ind_j]/strength_val[agents_ind_j]) * 
          dose_num_j * price[agents_ind_j] * (1 - discount(agents_ind_j)) + infusion_cost_j;
      }
      else{
        if (loading_dose[agents_ind_j] == 1){
          infusion_cost_j = 0;
        } 
        else{
          infusion_cost_j = infusion_cost[agents_ind_j] * dose_num[agents_ind_j];
        }
        tc += (ceil(weight * dose_amt[agents_ind_j]/strength_val[agents_ind_j]) * 
          dose_num_j * price[agents_ind_j] * (1 - discount(agents_ind_j)) + 
          infusion_cost_j) * cycle_length/12; 
      }
    }
  }
  return tc;
}

//// General management cost
// [[Rcpp::export]]
double sim_mgmt_cost1C(double yrlen, double cost){
  return cost * yrlen;
}

//// Hospitalization costs
// [[Rcpp::export]]
double sim_hosp_cost1C(double &haq, double yrlen, arma::rowvec hosp_days, 
                                   arma::rowvec cost_pday){
  double hosp_cost = 0.0;
  if (haq < 0.5){
      hosp_cost = hosp_days(0) * cost_pday(0) * yrlen;
  }
  else if (haq >= 0.5 && haq < 1){
      hosp_cost = hosp_days(1) * cost_pday(1) * yrlen;
  }
  else if (haq >= 1 && haq < 1.5){
      hosp_cost = hosp_days(2) * cost_pday(2) * yrlen;
  }
  else if (haq >= 1.5 & haq < 2){
      hosp_cost = hosp_days(3) * cost_pday(3) * yrlen;
  } 
  else if (haq >= 2 & haq < 2.5){
      hosp_cost = hosp_days(4) * cost_pday(4) * yrlen;
  }
  else if (haq >= 2.5){
      hosp_cost = hosp_days(5) * cost_pday(5) * yrlen;
  }
  return hosp_cost;
}

//// Serious Infection Cost
// [[Rcpp::export]]
double sim_si_cost1C(int si, double yrlen, double cost){
  double si_cost = 0.0;
  if (si == 1){
      si_cost = cost;
  }
  return si_cost;
}

// Productivity Loss
// [[Rcpp::export]]
double sim_prod_loss1C(double &haq, double yrlen, double beta){
  return haq * beta * yrlen;
}

// Wailoo 2006 HAQ to Utility Conversion
// [[Rcpp::export]]
double sim_utility_wailoo1C(double age, double disease_duration,
                           double haq0, int male, 
                           double prev_dmards, double haq, 
                           arma::rowvec b){
  double xb = b[0] + b[1] * age + b[2] * disease_duration + b[3] * haq0 +
    b[4] * male + b[5] * prev_dmards + b[6] * haq;
  return 1/(1 + exp(-xb));
}

//' @export
// [[Rcpp::export]]
std::vector<double> sim_utility_wailooC(std::vector<int> sim, std::vector<int> id,
                                        std::vector<double> age, double disease_duration,
                                        std::vector<double> haq0, std::vector<int> male, 
                                        std::vector<double> prev_dmards, std::vector<double> haq, 
                                        std::vector<double> b_int,
                                        std::vector<double> b_age, std::vector<double> b_disease_duration,
                                        std::vector<double> b_haq0, std::vector<double> b_male, 
                                        std::vector<double> b_prev_dmards, std::vector<double> b_haq
){
  int N = id.size();
  std::vector<double> logit_xb;
  logit_xb.reserve(N);
  double xb = 0.0;
  for (int i = 0; i < N; ++i){
    xb = b_int[sim[i]] + b_age[sim[i]] * age[i] + 
      b_disease_duration[sim[i]] * disease_duration + b_haq0[sim[i]] * haq0[id[i]] +
      b_male[sim[i]] * male[id[i]] + b_prev_dmards[sim[i]] * prev_dmards[id[i]] + 
      b_haq[sim[i]] * haq[i];
    logit_xb.push_back(1/(1 + exp(-xb)));
  }
  return logit_xb;
}

// Sample from Hernandez Alva (2013) Mixture Model
// Note the class 4 is the reference category in the paper 
//' @export
// [[Rcpp::export]]
double sim_utility_mixture1C(double haq,
                             double pain_mean, double haq_mean,
                             double pain_var, double haq_var, double painhaq_cor,
                             double age, int male,
                             arma::rowvec beta1, arma::rowvec beta2, 
                             arma::rowvec beta3, arma::rowvec beta4,
                             double alpha1, double alpha2,
                             double alpha3, double alpha4,
                             double alpha, 
                             double epsilon1_sd, double &epsilon2_sd,
                             double epsilon3_sd, double &epsilon4_sd,
                             double mu, arma::mat delta,
                             arma::rowvec w, arma::rowvec x){
  
  // update design matrices
  double pain = rbvcnormC(haq, pain_mean, haq_mean, pain_var, 
                           haq_var, painhaq_cor);
  x(0) = haq; x(1) = pow(haq, 2); x(2) = pain/100; 
  x(3) = age/10000; x(4) = pow(age/10000, 2);
  w(0) = 1.0; w(1) = haq; w(2) = pain/100; w(3) = pow(pain/100, 2);
  
  // latent class
  arma::rowvec latclass_prob = mlogit_probC(w, delta);
  int latclass = hesim::rcat1C(latclass_prob); 
  
  // utility conditional on latent class
  double ystar = 0.0;
  if (latclass == 0){
    double epsilon = R::rnorm(0, epsilon4_sd);
    ystar = alpha4 + alpha * male + dot(x, beta4) + mu + epsilon;
  }
  else if (latclass == 1){
    double epsilon = R::rnorm(0, epsilon1_sd);
    ystar = alpha1 + alpha * male + dot(x, beta1) + mu + epsilon;
  }
  else if (latclass == 2){
    double epsilon = R::rnorm(0, epsilon2_sd);
    ystar = alpha2 + alpha * male + dot(x, beta2) + mu + epsilon;
  }
  else {
    double epsilon = R::rnorm(0, epsilon3_sd);
    ystar = alpha3 + alpha * male + dot(x, beta3) + mu + epsilon;
  }
  if (ystar >= .883){
    ystar = 1.0;
  }
  return(ystar);
}

// Mixture Model HAQ to Utility Conversion
//' @export
// [[Rcpp::export]]
std::vector<double> sim_utility_mixtureC(std::vector<int> id, std::vector<int> sim, std::vector<double> haq,
                          double pain_mean, double haq_mean,
                          double pain_var, double haq_var, double painhaq_cor,
                          std::vector<double> age, std::vector<int> male,
                          arma::mat beta1, arma::mat beta2, 
                          arma::mat beta3, arma::mat beta4, 
                          std::vector<double> alpha1, std::vector<double> alpha2,
                          std::vector<double> alpha3, std::vector<double> alpha4,
                          std::vector<double> alpha,
                          std::vector<double> epsilon1_sd, std::vector<double> epsilon2_sd,
                          std::vector<double> epsilon3_sd, std::vector<double> epsilon4_sd,
                          std::vector<double> mu_sd, arma::cube delta){
  int N = haq.size();
  std::vector<double> util;
  util.reserve(N);
  arma::rowvec x(5);
  arma::rowvec w(4);
  int s = 0;
  double mu = R::rnorm(0, mu_sd[0]);
  for (int i = 0; i < N; ++i){
    s = sim[i];
    //sample for each new individual
    if (i > 0){
      if (sim[i] != sim[i-1] || id[i] != id[i-1]){
        mu = R::rnorm(0, mu_sd[s]);
      }
    }
    
    // sample for each individual at each time period
    util.push_back(sim_utility_mixture1C(haq[i],
                                         pain_mean, haq_mean, pain_var, haq_var, painhaq_cor,
                                         age[id[i]], male[id[i]],
                                         beta1.row(s), beta2.row(s),
                                         beta3.row(s), beta4.row(s),
                                         alpha1[s], alpha2[s], alpha3[s], alpha4[s], alpha[s],
                                        epsilon1_sd[s], epsilon2_sd[s], epsilon3_sd[s], epsilon4_sd[s],
                                        mu, delta.slice(s), w, x));
  }
  return util;
}

// Calculate means by individual for selected variables
class SimMeans {         
private:
  int sim; 
  int mod;
  int index;
  int n_indivs;
  int n_sims;
  int n_mods;
  double year;
  double discount_qalys;
  double discount_cost;
  std::map<std::string, std::vector<int> > id;
  std::map<std::string, std::vector<double> > varsums;
public:  
  SimMeans(int n_mods, int n_sims_, int n_indivs_, double r_qalys = .03, double r_cost = .03);
  std::map<std::string, std::vector<double> > get_varsums();
  std::map<std::string, std::vector<int> > get_id();
  void set_iterators(int m, int s, double y);
  void set_id();
  void increment_id(int m, int s, double y);
  void increment_varsums(double qalys); 
  std::map<std::string, std::vector<double> > calc_means();
};

SimMeans::SimMeans(int n_mods_, int n_sims_, int n_indivs_, double r_qalys, double r_cost){
  sim = 0;
  mod = 0;
  year = 0.0;
  n_indivs = n_indivs_;
  n_sims = n_sims_;
  n_mods = n_mods_;
  index = 0;
  discount_qalys = r_qalys;
  discount_cost = r_cost;
  int N = n_sims * n_mods;
  id = std::map<std::string, std::vector<int> > ();
  varsums = std::map<std::string, std::vector<double> >();
  id["sim"] = std::vector<int> (N);
  id["mod"] = std::vector<int> (N);
  varsums["qalys"] = std::vector<double> (N);
  varsums["dqalys"] = std::vector<double> (N);
}

std::map<std::string, std::vector<double> > SimMeans::get_varsums(){
  return varsums;
}

std::map<std::string, std::vector<int> > SimMeans::get_id(){
  return id;
}

void SimMeans::set_iterators(int m, int s, double y){
  mod = m;
  sim = s;
  year = y;
  index = mod * n_sims + sim;
}

void SimMeans::set_id(){
  id["mod"][index] = mod;
  id["sim"][index] = sim;
}

void SimMeans::increment_id(int m, int s, double y){
  set_iterators(m, s, y);
  set_id();
}

void SimMeans::increment_varsums(double qalys){
  varsums["qalys"][index] = varsums["qalys"][index] + qalys;
  varsums["dqalys"][index] = varsums["dqalys"][index] + qalys * discount_factor(year, discount_qalys);
}

std::map<std::string, std::vector<double> > SimMeans::calc_means(){
  std::map<std::string, std::vector<double> > varmeans(varsums);
  std::map<std::string, std::vector<double> >::iterator it;
  for (it = varmeans.begin(); it != varmeans.end(); ++it){
    int J = it->second.size();
    for (int j = 0; j < J; ++j){
      it->second[j] = it->second[j]/n_indivs;
    }
  }
  return varmeans;
}

RCPP_MODULE(mod_SimMeans) {
    
    class_<SimMeans>("SimMeans")
    .constructor<int, int, int, double, double>()
    .method("get_id", &SimMeans::get_id)
    .method("get_varsums", &SimMeans::get_varsums)
    .method("set_iterators", &SimMeans::set_iterators)
    .method("set_id", &SimMeans::set_id)
    .method("increment_id", &SimMeans::increment_id)
    .method("increment_varsums", &SimMeans::increment_varsums)
    .method("calc_means", &SimMeans::calc_means)
    ;
}

// Calculate means by time period for selected variables
class TimeMeans {         
private:
  int sim, mod;
  double month, cycle_length;
  int index;
  int n_indivs, n_sims, n_mods, n_cycles;
  std::vector<int> alive;
  std::map<std::string, std::vector<int> > id;
  std::map<std::string, std::vector<double> > varsums;
public:  
  TimeMeans(int n_mods_, int n_sims_, int n_indivs_, int ncycles_, double cycle_length_);
  std::map<std::string, std::vector<double> > get_varsums();
  std::map<std::string, std::vector<int> > get_id();
  int get_index();
  std::vector<int> get_alive();
  void set_iterators(int m, int s, double month_);
  void set_id();
  void increment_id(int m, int s, double month_);
  void increment_alive();
  void increment_varsums(double qalys, double haq); 
  std::map<std::string, std::vector<double> > calc_means();
};

TimeMeans::TimeMeans(int n_mods_, int n_sims_, int n_indivs_, int n_cycles_, 
                     double cycle_length_){
  sim = 0;
  mod = 0;
  month = 0.0;
  cycle_length = cycle_length_;
  n_indivs = n_indivs_;
  n_sims = n_sims_;
  n_mods = n_mods_;
  n_cycles = n_cycles_;
  index = 0;
  int N = n_sims * n_mods * n_cycles;
  alive = std::vector<int> (N);
  std::fill (alive.begin(), alive.end(), 0);
  id["sim"] = std::vector<int> (N);
  id["mod"] = std::vector<int> (N);
  id["month"] = std::vector<int> (N);
  varsums["qalys"] = std::vector<double> (N);
  varsums["haq"] = std::vector<double> (N);
}

std::map<std::string, std::vector<double> > TimeMeans::get_varsums(){
  return varsums;
}

std::map<std::string, std::vector<int> > TimeMeans::get_id(){
  return id;
}

int TimeMeans::get_index(){
  return index;
}

std::vector<int> TimeMeans::get_alive(){
  return alive;
}

void TimeMeans::set_iterators(int m, int s, double month_){
  mod = m;
  sim = s;
  month = month_;
  int cycle = (int) month/cycle_length - 1;
  index = mod * n_sims * n_cycles +  sim * n_cycles + cycle;
}

void TimeMeans::set_id(){
  id["mod"][index] = mod;
  id["sim"][index] = sim;
  id["month"][index] = month;
}

void TimeMeans::increment_id(int m, int s, double month_){
  set_iterators(m, s, month_);
  set_id();
}

void TimeMeans::increment_alive(){
  alive[index] += 1;
}

void TimeMeans::increment_varsums(double qalys, double haq){
  varsums["qalys"][index] = varsums["qalys"][index] + qalys;
  varsums["haq"][index] = varsums["haq"][index] + haq;
}

std::map<std::string, std::vector<double> > TimeMeans::calc_means(){
  std::map<std::string, std::vector<double> > varmeans(varsums);
  std::map<std::string, std::vector<double> >::iterator it;
  for (it = varmeans.begin(); it != varmeans.end(); ++it){
    int J = it->second.size();
    for (int j = 0; j < J; ++j){
      if (alive[j] > 0){
        it->second[j] = it->second[j]/alive[j];
      }
      else{
        it->second[j] = NA_REAL;
      }
    }
  }
  return varmeans;
}

RCPP_MODULE(mod_TimeMeans) {
  
  class_<TimeMeans>("TimeMeans")
  .constructor<int, int, int, int, double>()
  .method("get_varsums", &TimeMeans::get_varsums)
  .method("get_id", &TimeMeans::get_id)
  .method("get_index", &TimeMeans::get_index)
  .method("get_alive", &TimeMeans::get_alive)
  .method("increment_id", &TimeMeans::increment_id)
  .method("increment_alive", &TimeMeans::increment_alive)
  .method("increment_varsums", &TimeMeans::increment_varsums)
  .method("calc_means", &TimeMeans::calc_means)
  ;
}

// Output at time-period 0 for each treatment 
class Out0 {         
private:
  int n_indivs, n_sims, n_mods, n_tx;
  std::vector<int> sim_vec, mod_vec, id_vec, tx_vec;
  std::vector<int> acr_vec, eular_vec;
  std::vector<double> ttd_vec, ttsi_vec;
public:  
  Out0(int n_mods_, int n_sims_, int n_indivs_, int n_tx_);
  std::vector<int> get_sim();
  std::vector<int> get_mod();
  std::vector<int> get_id();
  std::vector<int> get_tx();
  std::vector<int> get_acr();
  std::vector<int> get_eular();
  std::vector<double> get_ttd();
  std::vector<double> get_ttsi();
  void push_back(int cycle, int mod, int sim, int id, int tx, int acr, int eular, 
                 double ttd, double ttsi);
};

Out0::Out0(int n_mods_, int n_sims_, int n_indivs_, int n_tx_){
  n_indivs = n_indivs_;
  n_sims = n_sims_;
  n_mods = n_mods_;
  n_tx = n_tx_;
  int N = n_sims * n_mods * n_indivs * n_tx;
  sim_vec.reserve(N);
  mod_vec.reserve(N);
  id_vec.reserve(N);
  tx_vec.reserve(N);
  acr_vec.reserve(N);
  eular_vec.reserve(N);
  ttd_vec.reserve(N);
  ttsi_vec.reserve(N);
}

std::vector<int> Out0::get_sim(){
  return sim_vec;
}

std::vector<int> Out0::get_mod(){
  return mod_vec;
}

std::vector<int> Out0::get_id(){
  return id_vec;
}

std::vector<int> Out0::get_tx(){
  return tx_vec;
}

std::vector<int> Out0::get_acr(){
  return acr_vec;
}

std::vector<int> Out0::get_eular(){
  return eular_vec;
}

std::vector<double> Out0::get_ttd(){
  return ttd_vec;
}

std::vector<double> Out0::get_ttsi(){
  return ttsi_vec;
}

void Out0::push_back(int cycle, int mod, int sim, int id, int tx, int acr, int eular, 
                     double ttd, double ttsi){
  if (cycle == 0){
    mod_vec.push_back(mod);
    sim_vec.push_back(sim);
    id_vec.push_back(id);
    tx_vec.push_back(tx);
    acr_vec.push_back(acr);
    eular_vec.push_back(eular);
    ttd_vec.push_back(ttd);
    ttsi_vec.push_back(ttsi);
  }
}

RCPP_MODULE(mod_Out0) {
  
  class_<Out0>("Out0")
  .constructor<int, int, int, int>()
  .method("get_sim", &Out0::get_sim)
  .method("get_mod", &Out0::get_mod)
  .method("get_id", &Out0::get_id)
  .method("get_tx", &Out0::get_tx)
  .method("get_acr", &Out0::get_acr)
  .method("get_eular", &Out0::get_eular)
  .method("get_ttd", &Out0::get_ttd)
  .method("get_ttsi", &Out0::get_ttsi)
  .method("push_back", &Out0::push_back)
  ;
}

// Calculate vector of QALYs given simulation output
//' @export
// [[Rcpp::export]]
std::vector<double> qalysC(std::vector<double> &utility, std::vector<double> &yrlen,
                           std::vector<int> &sim, std::vector<int> &si, std::vector<double> &si_ul){
  int N = utility.size();
  std::vector<double> qalys_vec;
  qalys_vec.reserve(N);
  for (int i = 0; i < N; ++i){
    qalys_vec.push_back(yrlen[i] * (utility[i] - si[i] * si_ul[sim[i]]/12));
  }
  return qalys_vec;
}

// Simulate HAQ score
// [[Rcpp::export]] 
List sim_iviRA_C(arma::mat arm_inds, CharacterMatrix model_structures_mat,
             std::vector<double> haq0, std::vector<double> das28_0,
             std::vector<double> sdai0, std::vector<double> cdai0,
             std::vector<double> age0, std::vector<int> male,
             std::vector<int> prev_dmards,
             arma::cube nma_acr1, arma::cube nma_acr2, arma::mat nma_dhaq1, arma::mat nma_dhaq2,
             arma::mat nma_das28_1, arma::mat nma_das28_2,
             arma::cube acr2eular,arma::mat acr2haq, arma::mat eular2haq, 
             arma::mat acr2das28, arma::mat acr2sdai, arma::mat acr2cdai,
             arma::mat tswitch_da,
             arma::mat haq_lprog_therapy, arma::mat haq_lprog_age,
             arma::cube haq_lcgm_delta, arma::cube haq_lcgm_beta,
             std::vector<double> rebound_factor,
             arma::mat lifetable_male, arma::mat lifetable_female, 
             arma::mat x_mort, arma::mat logor_mort, 
             arma::mat x_ttd_all, arma::mat x_ttd_da, arma::mat x_ttd_eular, 
             Rcpp::List ttd_all_list, Rcpp::List ttd_da_list, Rcpp::List ttd_eular_mod_list, Rcpp::List ttd_eular_good_list,
             int cdmards, int nbt, 
             arma::mat si_loc, arma::mat si_anc1, arma::mat si_anc2, std::string si_dist, 
             arma::mat haqdelta_loghr, int max_months, 
             arma::mat hosp_days, arma::mat cost_pday, std::vector<double> mgmt_cost, 
             std::vector<double> si_cost, std::vector<double> prod_loss, 
             Rcpp::List tc_list, std::vector<double> weight, 
             arma::mat coefs_wailoo, Rcpp::List pars_util_mix, std::vector<double> si_ul,
             Rcpp::List discount_rate, std::string output){
  
  // Type conversions
  //// Model structures
  std::vector<std::vector<std::string> > model_structures = charmat2stdvec(model_structures_mat);
  
  //// Treatment costs
  arma::cube tc_agents_ind = as<arma::cube>(tc_list["agents"]);
  arma::mat tc_discount = as<arma::mat>(tc_list["discount"]);
  Rcpp::DataFrame tc = Rcpp::as<Rcpp::DataFrame>(tc_list["cost"]);
  std::vector<std::string> tx_name = as<std::vector<std::string> >(tc["sname"]);
  std::vector<double> init_dose_val = as<std::vector<double> >(tc["init_dose_val"]);
  std::vector<double> ann_dose_val = as<std::vector<double> >(tc["ann_dose_val"]);
  std::vector<double> strength_val = as<std::vector<double> >(tc["strength_val"]);
  std::vector<double> init_num_doses = as<std::vector<double> >(tc["init_num_doses"]);
  std::vector<double> ann_num_doses = as<std::vector<double> >(tc["ann_num_doses"]);
  std::vector<double> price = as<std::vector<double> >(tc["price_per_unit"]);
  std::vector<double> infusion_cost = as<std::vector<double> >(tc["infusion_cost"]);
  std::vector<int> loading_dose = as<std::vector<int> >(tc["loading_dose"]);
  std::vector<int> weight_based = as<std::vector<int> >(tc["weight_based"]);
  
  //// Utility mixture model
  arma::mat utilmix_beta1 = as<arma::mat> (pars_util_mix["beta1"]);
  arma::mat utilmix_beta2 = as<arma::mat> (pars_util_mix["beta2"]);
  arma::mat utilmix_beta3 = as<arma::mat> (pars_util_mix["beta3"]);
  arma::mat utilmix_beta4 = as<arma::mat> (pars_util_mix["beta4"]);
  std::vector<double> utilmix_alpha1 = as<std::vector<double> > (pars_util_mix["alpha1"]);
  std::vector<double> utilmix_alpha2 = as<std::vector<double> > (pars_util_mix["alpha2"]);
  std::vector<double> utilmix_alpha3 = as<std::vector<double> > (pars_util_mix["alpha3"]);
  std::vector<double> utilmix_alpha4 = as<std::vector<double> > (pars_util_mix["alpha4"]);
  std::vector<double> utilmix_alpha = as<std::vector<double> > (pars_util_mix["alpha"]);
  std::vector<double> utilmix_epsilon1_sd = as<std::vector<double> > (pars_util_mix["epsilon1"]);
  std::vector<double> utilmix_epsilon2_sd = as<std::vector<double> > (pars_util_mix["epsilon2"]);
  std::vector<double> utilmix_epsilon3_sd = as<std::vector<double> > (pars_util_mix["epsilon3"]);
  std::vector<double> utilmix_epsilon4_sd = as<std::vector<double> > (pars_util_mix["epsilon4"]);
  std::vector<double> utilmix_mu_sd = as<std::vector<double> > (pars_util_mix["mu"]);
  arma::cube utilmix_delta = as<arma::cube> (pars_util_mix["delta"]);
  double pain_mean = as<double>(as<Rcpp::List>(pars_util_mix["pain"])["pain.mean"]);
  double haq_mean = as<double>(as<Rcpp::List>(pars_util_mix["pain"])["haq.mean"]);
  double pain_var = as<double>(as<Rcpp::List>(pars_util_mix["pain"])["pain.var"]);
  double haq_var = as<double>(as<Rcpp::List>(pars_util_mix["pain"])["haq.var"]);
  double painhaq_cor = as<double>(as<Rcpp::List>(pars_util_mix["pain"])["painhaq.cor"]);
  
  //// Time to treatment discontinuation
  TTDPars ttd_all = get_ttd_pars(ttd_all_list);
  TTDPars ttd_da = get_ttd_pars(ttd_da_list);
  TTDPars ttd_eular_mod = get_ttd_pars(ttd_eular_mod_list);
  TTDPars ttd_eular_good = get_ttd_pars(ttd_eular_good_list);
  
  //// Discount rate
  double discount_qalys = as<double>(discount_rate["qalys"]);
  double discount_cost = as<double>(discount_rate["cost"]);
  
  // Constants
  const double cycle_length = 6;
  const int maxt = 100 * (12/cycle_length);
  const int max_cycles = (int) max_months/cycle_length;
  const int treat_gap = 0;
  const int n_mods = model_structures.size();
  const int n_treatments = arm_inds.n_cols;
  const int n_sims = logor_mort.n_rows;
  const int n_pat = x_mort.n_rows;
  const double disease_duration = 18.65;
  
  // Declarations
  int counter = 0;
  double cage_i = 0.0; // continuous age
  int age_i = 0;
  int arm_row = 0;
  arma::rowvec arm_ind_i = arm_inds.row(0);
  int arm_ind_ij = 0;
  Cost sim_cost;
  ModelStructure mod_struct;
  double utility = 0.0;
  arma::rowvec utilmix_x(5);
  arma::rowvec utilmix_w(4);
  double qalys = 0.0;
  SimMeans sim_means(n_mods, n_sims, n_pat, discount_qalys, discount_cost);
  TimeMeans time_means(n_mods, n_sims, n_pat, std::min(maxt, max_cycles),
                       cycle_length);
  Out0 out0(n_mods, n_sims, n_pat, n_treatments);
  
  // Information to store
  std::vector<int> model;
  std::vector<int> sim;
  std::vector<int> id;
  std::vector<double> month_vec;
  std::vector<int> tx_vec;
  std::vector<int> tx_seq;
  std::vector<int> tx_cycle;
  std::vector<int> death_vec;
  std::vector<double> age_vec;
  std::vector<double> ttd_vec;
  std::vector<int> acr_vec;
  std::vector<int> eular_vec;
  std::vector<double> das28_vec;
  std::vector<double> sdai_vec;
  std::vector<double> cdai_vec;
  std::vector<double> haq_vec;
  std::vector<double> ttsi_vec;
  std::vector<int> si_vec;
  std::vector<double> yrlen_vec;
  std::vector<double> tx_cost_vec;
  std::vector<double> hosp_cost_vec;
  std::vector<double> mgmt_cost_vec;
  std::vector<double> si_cost_vec;
  std::vector<double> prod_loss_vec;
  std::vector<double> utility_vec;
  std::vector<double> qalys_vec;
  
  // Loop over models
  for (int m = 0; m < n_mods; ++m){
    mod_struct = update_model_structure(mod_struct, model_structures[m]);
    
    // "OUTER" LOOP
    for (int s = 0; s < n_sims; ++s){
      
      // "INNER" LOOP
      // Loop over Individuals
      for (int i = 0; i < n_pat; ++i){
        double month = 6;
        cage_i = age0[i];
        age_i = int(cage_i);
        double haq = haq0[i];
        if (arm_inds.n_rows > 1){
          arm_row = i;
          arm_ind_i = arm_inds.row(i);
        }
        int t_cdmards = 0;
        Tx_ISwitch sim_s_t1;
        sim_s_t1.das28 = das28_0[i];
        sim_s_t1.sdai = sdai0[i];
        sim_s_t1.cdai = cdai0[i];
        arma::rowvec x_ttd_da_i;
        double utilmix_mu = R::rnorm(0, utilmix_mu_sd[s]);
        
        // Loop over treatments
        for (int j = 0; j < n_treatments + 1; ++j){
          int si = 0;
          
          // Which treatment is being used and how long has it been used for
          if (j == n_treatments - 1){
             t_cdmards = 0; // counter for cdmards and nbt which does not reset after cdmards ends and nbt begins
          }
          
          if (j == n_treatments){
            arm_ind_ij = nbt;
            if (arm_ind_i(j-1) == cdmards){
                t_cdmards = tx_cycle[counter - 1] + 1; // note 0 otherwise
            } 
          }
          else{
            arm_ind_ij = arm_ind_i(j);
          }
          
          // H1-H3: simulate change in HAQ during initial treatment phase
          TxIHaq sim_h_t1 = sim_tx_ihaq(mod_struct.tx_ihaq, j, arm_ind_ij, nbt,
                                          nma_acr1.slice(arm_ind_ij).row(s), 
                                          nma_acr2.slice(arm_ind_ij).row(s), 
                                          nma_dhaq1.col(arm_ind_ij)(s),
                                          nma_dhaq2.col(arm_ind_ij)(s),
                                          acr2eular.slice(s), acr2haq.row(s), 
                                          eular2haq.row(s));
          
          // S1-S6: simulate treatment switching during first 6 months
          sim_s_t1 = sim_tx_iswitch(mod_struct.tx_iswitch, j, arm_ind_ij, nbt,
                                                    sim_h_t1.acr, sim_h_t1.eular,
                                                    sim_s_t1.das28, sim_s_t1.sdai, sim_s_t1.cdai,
                                                    acr2das28.row(s), acr2sdai.row(s), acr2cdai.row(s),
                                                    nma_das28_1.row(s), nma_das28_2.row(s), 
                                                    tswitch_da.row(s));
          
          // Time to treatment discontinuation and serious infections
          double ttsi_j = (-cycle_length/12 + rsurvC(as_scalar(si_loc.row(s).col(arm_ind_ij)), 
                                                as_scalar(si_anc1.row(s).col(arm_ind_ij)),
                                                si_dist,as_scalar(si_anc2.row(s).col(arm_ind_ij))
                                                )) * (12/cycle_length);
          double ttd_j = 0;
          if (mod_struct.tx_iswitch == "acr-eular-switch"){
              ttd_j = sim_ttd_eular(x_ttd_eular.row(i), 
                                    ttd_eular_mod.loc[mod_struct.ttd_dist].row(s), ttd_eular_mod.anc1[mod_struct.ttd_dist](s), 
                                    ttd_eular_good.loc[mod_struct.ttd_dist].row(s), ttd_eular_good.anc1[mod_struct.ttd_dist](s),
                                    sim_h_t1.eular, mod_struct.ttd_dist, cycle_length, ttsi_j,
                                    ttd_eular_mod.anc2[mod_struct.ttd_dist](s), ttd_eular_good.anc2[mod_struct.ttd_dist](s));
          }
          else if (mod_struct.tx_iswitch == "acr-switch"){
            ttd_j = sim_ttd(x_ttd_all.row(i), ttd_all.loc[mod_struct.ttd_dist].row(s), ttd_all.anc1[mod_struct.ttd_dist](s),
                                sim_s_t1.tswitch, mod_struct.ttd_dist, cycle_length, ttsi_j,
                                ttd_all.anc2[mod_struct.ttd_dist](s));
          }
          else {
            x_ttd_da_i = update_x_ttd_da(x_ttd_da.row(i), sim_s_t1.da_cat);
            ttd_j = sim_ttd(x_ttd_da.row(i), ttd_da.loc[mod_struct.ttd_dist].row(s), ttd_da.anc1[mod_struct.ttd_dist](s),
                            sim_s_t1.tswitch, mod_struct.ttd_dist, cycle_length, ttsi_j,
                            ttd_da.anc2[mod_struct.ttd_dist](s));
          }
          
          // Loop over time
          for (int t = 0; t < maxt; ++t){
            
            // Stop j loop if cycle is larger than treatment duration
            if (j < n_treatments && t > 0 && t >= ttd_j + 1 + treat_gap){
              cage_i = cage_i - cycle_length/12 + 0.5;
              month = month - cycle_length + 6;
              age_i = int(cage_i);
              break;
            }
            
            // Draw death indicator
            double death = sample_deathC(age_i, male[i], lifetable_male, lifetable_female, 
                                         x_mort.row(i), logor_mort.row(s), haq0[i], haq,
                                        cycle_length, month, haqdelta_loghr.row(s));
            
            // Did serious infection occur?
            if (ttsi_j < ttd_j && t > ttd_j && t < ttd_j + 1 && arm_ind_ij != nbt) {
              si = 1;
            } 
            else if (ttsi_j < 0 && arm_ind_ij != nbt) {
             si = 1;
            }
            
            // Update HAQ score
            if (t == 0){ // initial haq change
              if (ttd_j == 0){ // immediate rebound for patients switching treatment during the initial period
                  // haq = haq + sim_h_t1.dhaq - sim_h_t1.dhaq * rebound_factor[s]; 
                  haq = haq;
              } else{
                haq = haq + sim_h_t1.dhaq; 
              }
            }
            else if (t > 0 && t <= ttd_j){
              if(mod_struct.cdmards_haq_model == "lcgm" && (arm_ind_ij == nbt || arm_ind_ij == cdmards)){
                  haq = haq + sim_dhaq_lcgm1C(2.5 + t_cdmards * cycle_length/12, cycle_length, cage_i, 1 - male[i],
                                             das28_0[i], haq_lcgm_delta.slice(s), haq_lcgm_beta.slice(s));
              } else{
                update_haq_t(haq, as_scalar(haq_lprog_therapy.row(s).col(arm_ind_ij)),
                             haq_lprog_age.row(s), cage_i, cycle_length);
              }
            }
            else if (t > ttd_j && t < ttd_j + 1 && j < n_treatments){ // rebound
              haq = haq - sim_h_t1.dhaq * rebound_factor[s];
            }
            else {
              update_haq_t(haq, as_scalar(haq_lprog_therapy.row(s).col(nbt)),
                            haq_lprog_age.row(s), cage_i, cycle_length); // consider making j nbt here as patient is no longer on treatment
            }
            if (haq < 0){
              haq = 0.0;
            }
            if (haq > 3){
              haq = 3.0;
            }
            
            // Simulate costs
            sim_cost.tx = sim_tx_cost1C(t, tc_agents_ind.slice(arm_row).row(j), tx_name,
                                        init_dose_val, ann_dose_val, strength_val,
                                        init_num_doses, ann_num_doses, price,
                                        infusion_cost, loading_dose, 
                                        weight_based, weight[i], cycle_length, 
                                        tc_discount.row(s));
            sim_cost.hosp = sim_hosp_cost1C(haq, cycle_length/12, hosp_days.row(s), cost_pday.row(s));
            sim_cost.mgmt = sim_mgmt_cost1C(cycle_length/12, mgmt_cost[s]);
            sim_cost.prod = sim_prod_loss1C(haq, cycle_length/12, prod_loss[s]);
            sim_cost.si = sim_si_cost1C(si, cycle_length/12, si_cost[s]);
            
            // Simulate utility + QALYs
            if (mod_struct.utility_model == "wailoo"){
              utility = sim_utility_wailoo1C(cage_i, disease_duration, haq0[i],  male[i],
                                            prev_dmards[i], haq, coefs_wailoo.row(s));
            }
            else if (mod_struct.utility_model == "mixture"){
                 utility = sim_utility_mixture1C(haq, pain_mean, haq_mean,
                                                 pain_var, haq_var, painhaq_cor,
                                                 cage_i, male[i], 
                                                 utilmix_beta1.row(s), utilmix_beta2.row(s),
                                                 utilmix_beta3.row(s), utilmix_beta4.row(s),
                                                 utilmix_alpha1[s], utilmix_alpha2[s],
                                                 utilmix_alpha3[s], utilmix_alpha4[s],
                                                 utilmix_alpha[s], 
                                                 utilmix_epsilon1_sd[s], utilmix_epsilon2_sd[s],
                                                 utilmix_epsilon3_sd[s], utilmix_epsilon4_sd[s],
                                                 utilmix_mu, utilmix_delta.slice(s),
                                                 utilmix_w, utilmix_x);
            }
            qalys = cycle_length/12 * (utility - si * si_ul[s]);
            
            if (output == "summary"){
              sim_means.increment_id(m, s, month/12);
              sim_means.increment_varsums(qalys);
              time_means.increment_id(m, s, month); //note: requires assumption of constant model cycles!!!
              time_means.increment_alive();
              time_means.increment_varsums(qalys, haq);
              out0.push_back(t, m, s, i, j, sim_h_t1.acr, sim_h_t1.eular,
                             ttd_j, ttsi_j);
            }
            
            // Fill vectors
            model.push_back(m);
            sim.push_back(s);
            id.push_back(i);
            month_vec.push_back(month);
            tx_vec.push_back(arm_ind_ij);
            tx_seq.push_back(j);
            tx_cycle.push_back(t);
            death_vec.push_back(death);
            age_vec.push_back(cage_i);
            if (arm_ind_ij != nbt){
              ttd_vec.push_back(ttd_j - t);
            }
            else {
              ttd_vec.push_back(NA_REAL);
            }
            
            acr_vec.push_back(sim_h_t1.acr);
            eular_vec.push_back(sim_h_t1.eular);
            das28_vec.push_back(sim_s_t1.das28);
            sdai_vec.push_back(sim_s_t1.sdai);
            cdai_vec.push_back(sim_s_t1.cdai);
            haq_vec.push_back(haq);
            if (arm_ind_ij != nbt){
              ttsi_vec.push_back(ttsi_j - t);
            }
            else{
              ttsi_vec.push_back(NA_REAL);
            }
            si_vec.push_back(si);
            if (t == 0) {
              yrlen_vec.push_back(0.5);
            } 
            else {
              yrlen_vec.push_back(cycle_length/12);
            }
            tx_cost_vec.push_back(sim_cost.tx);
            hosp_cost_vec.push_back(sim_cost.hosp);
            mgmt_cost_vec.push_back(sim_cost.mgmt);
            prod_loss_vec.push_back(sim_cost.prod);
            si_cost_vec.push_back(sim_cost.si);
            utility_vec.push_back(utility);
            qalys_vec.push_back(qalys);
            
            // Iterate and update age + time
            cage_i = cage_i + cycle_length/12;
            month = month + cycle_length;
            age_i = int(cage_i);
            t_cdmards = t_cdmards + 1;
            ++counter;
            
            // Stop i loop for death
            if (death == 1 || max_months < month){
              goto next_i;
            }
          }
        } next_i:;
      }
    }
  }
  // Return results
  if (output == "data"){
    Rcpp::DataFrame sim1 = Rcpp::DataFrame::create(
      _["model"] = model,
      _["sim"] = sim, 
      _["id"] = id,
      _["month"] = month_vec,
      _["tx"] = tx_vec,
      _["tx_seq"] = tx_seq,
      _["tx_cycle"] = tx_cycle,
      _["death"] = death_vec,
      _["age"] = age_vec,
      _["ttd"] = ttd_vec,
      _["acr"] = acr_vec,
      _["eular"] = eular_vec,
      _["das28"] = das28_vec,
      _["sdai"] = sdai_vec,
      _["cdai"] = cdai_vec,
      _["haq"] = haq_vec,
      _["ttsi"] = ttsi_vec,
      _["si"] = si_vec
    );
    Rcpp::DataFrame sim2 = Rcpp::DataFrame::create(
      _["yrlen"] = yrlen_vec,
      _["tx_cost"] = tx_cost_vec,
      _["hosp_cost"] = hosp_cost_vec,
      _["mgmt_cost"] = mgmt_cost_vec,
      _["si_cost"] = si_cost_vec,
      _["prod_loss"] = prod_loss_vec,
      _["utility"] = utility_vec,
      _["qalys"] = qalys_vec
    );
    return Rcpp::List::create(sim1, sim2);
  } 
  else{
    std::map<std::string, std::vector<double> > sim_means_out = sim_means.calc_means();
    std::map<std::string, std::vector<int> > sim_means_id = sim_means.get_id();
    Rcpp::DataFrame sim_means_df = Rcpp::DataFrame::create(
      _["model"] = sim_means_id["mod"],
      _["sim"] = sim_means_id["sim"],
      _["qalys"] = sim_means_out["qalys"],
      _["dqalys"] = sim_means_out["dqalys"] 
    );
    std::map<std::string, std::vector<double> > time_means_out = time_means.calc_means();
    std::map<std::string, std::vector<int> > time_means_id = time_means.get_id();
    Rcpp::DataFrame time_means_df = Rcpp::DataFrame::create(
      _["model"] = time_means_id["mod"],
      _["sim"] = time_means_id["sim"],
      _["month"] = time_means_id["month"],
      _["alive"] = time_means.get_alive(),
      _["qalys"] = time_means_out["qalys"],
      _["haq"] = time_means_out["haq"] 
    );
    Rcpp::DataFrame out0_df = Rcpp::DataFrame::create(
      _["model"] = out0.get_mod(),
      _["sim"] = out0.get_sim(),
      _["id"] = out0.get_id(),
      _["tx"] = out0.get_tx(),
      _["acr"] = out0.get_acr(),
      _["eular"] = out0.get_eular(),
      _["ttd"] = out0.get_ttd(),
      _["ttsi"] = out0.get_ttsi()
    );
    return Rcpp::List::create(
      _["means"] = sim_means_df,
      _["time.means"] = time_means_df,
      _["out0"] = out0_df
      );
  }
}




