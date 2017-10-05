// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(hesim)]]

#include <RcppArmadillo.h>
#include "rand.h"
#include "mort.h"
#include "logit-prob.h"
#include "utils.h"
#include "ips.h"
#include <hesim.h>
using namespace Rcpp;

/*****************
* Model structure
*****************/
struct ModelStructure {
  std::string tx_ihaq;
  std::string tx_iswitch;
  std::string cdmards_haq_model;
  std::string ttd_cause;
  std::string ttd_dist;
  std::string utility_model;
};

ModelStructure update_model_structure(ModelStructure m, std::vector<std::string> x){
  m.tx_ihaq = x[0];
  m.tx_iswitch = x[1];
  m.cdmards_haq_model = x[2];
  m.ttd_cause = x[3];
  m.ttd_dist = x[4];
  m.utility_model = x[5];
  return m;
}

/************************************
* Update HAQ with linear progression 
************************************/
void update_haq_t1(double &haq, double haq_change){
  haq = haq + haq_change;
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

/**************
* NMA ACR class
**************/
void nmaACR::set(std::string hist_, double rr_, double A_, double z2_, double z3_,
                 arma::rowvec d_beta_, arma::rowvec x_, int line){
  hist = hist_;
  if (line == 0 && hist == "naive"){
    rr = 1;
  } 
  else{
    rr = rr_;
  }
  A = A_;
  z2 = z2_;
  z3 = z3_;
  d_beta = d_beta_;
  x = x_;
}

arma::rowvec nmaACR::nma_acrprob(){
  double d = arma::dot(d_beta, x);
  
  // probability less than category
  double pl_70 = R::pnorm5(A + z3 + d, 0, 1, 1, 0); // less than ACR 70
  double pl_50 = R::pnorm5(A + z2 + d, 0, 1, 1, 0); // less than ACR 50
  double pl_20 = R::pnorm5(A + d, 0, 1, 1, 0); // less than ACR 20
  
  // probability in overlapping categories
  double po_g70 =  rr * (1 - pl_70); // greater than ACR 70
  double po_g50 = rr * (1 - pl_50); // greater than ACR 50
  double po_g20 = rr * (1 - pl_20); // greater than ACR 20
  double po_l20 = 1 - po_g20; // less than ACR 20
  
  // probability in mutually exclusive categories
  arma::rowvec p(4);
  p(0) = po_l20;  // less than ACR 20
  p(1) = po_g20 - po_g50; // ACR 20-50
  p(2) = po_g50 - po_g70; // ACR 50-70
  p(3) = po_g70; // greater than ACR 70
  
  // return
  return p;
}

double nmaACR::sim_acr(){
  arma::rowvec p = nma_acrprob();
  double acr = hesim::rcat1C(p);
  return acr;
}

/****************************
* NMA class for linear model
****************************/
void nmaLM::set(std::string hist_, double rr_, double A_,
                 arma::rowvec d_beta_, arma::rowvec x_, int line){
  hist = hist_;
  if (line == 0 && hist == "naive"){
    rr = 1;
  } 
  else{
    rr = rr_;
  }
  A = A_;
  d_beta = d_beta_;
  x = x_;
}

double nmaLM::sim_dy(){
  double dy = A + arma::dot(d_beta, x);
  if (dy <= 0){
    return rr * dy;
  } 
  else{
    return (1 + 1 - rr) * dy;
  }
}

/****************************************************************
* TxIHAQ class:
* Effect of treatment on HAQ during the initial treatment phase
****************************************************************/
void TxIHaq::sim(std::string model, int line, int therapy, int nbt,
                    nmaACR nma_acr, nmaLM nma_dhaq,
                     arma::mat acr2eular, arma::rowvec acr2haq, arma::rowvec eular2haq){
  
  // ACR response
  if (model == "acr-haq" || model == "acr-eular-haq"){
      if (therapy != nbt){
          acr = nma_acr.sim_acr();
      }
      else{
          acr = 0;
      }
  }
  
  // EULAR response
  if (model == "acr-eular-haq"){
      if (therapy != nbt){
        eular = hesim::rcat1C(acr2eular.row(acr));
      }
      else{
        eular = 0;
      }
  }
  
  // HAQ
  if (model == "acr-haq" && therapy != nbt){
        dhaq = acr2haq(acr);
  } 
  else if (model == "acr-eular-haq" && therapy != nbt){
      dhaq = eular2haq(eular);
  } 
  else if (model == "haq" && therapy != nbt) {
      dhaq = nma_dhaq.sim_dy();
  }
  else if (therapy == nbt){
      dhaq = 0.0;
  }
}

/**********************************************************************
* TxISwitch class:
* Effect of treatment on switching treatment during the treatment phase
**********************************************************************/
int TxISwitch::get_das28_cat(){
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

int TxISwitch::get_sdai_cat(){
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

int TxISwitch::get_cdai_cat(){
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

double TxISwitch::get_da_new(double da_old, double da_change, double lower, double upper){
  double da_new = da_old + da_change;
  if (da_new < lower){
    da_new = 0;
  } else if (da_new >upper){
    da_new = upper;
  }
  return da_new;
}

void TxISwitch::sim(std::string model, int line, int therapy, int nbt,
                  int acr, int eular,
                  arma::rowvec acr2das28, arma::rowvec acr2sdai,arma::rowvec acr2cdai,
                  nmaLM nma_dhaq, arma::rowvec p){
  
  tswitch = 0;
  double da_change = 0.0;
  
  // Switching rules
  if (model == "acr-switch"){ // Treatment -> ACR -> Switch
    if (acr == 0){
      tswitch = 1;
    }
  }
  else if (model == "acr-das28-switch" ||   // Treatment -> ACR -> DA -> Switch
           model == "acr-sdai-switch" ||
           model == "acr-cdai-switch") {
      if (model == "acr-das28-switch"){
        da_change = acr2das28(acr);
        das28 = get_da_new(das28, da_change, 0, 9.4);
        da_cat = get_das28_cat();
      }
      else if (model == "acr-sdai-switch"){
        da_change = acr2sdai(acr);
        sdai = get_da_new(sdai, da_change, 0, 86);
        da_cat = get_sdai_cat();
      }
      else if (model == "acr-cdai-switch"){
        da_change = acr2cdai(acr);
        cdai = get_da_new(cdai, da_change, 0, 76);
        da_cat = get_cdai_cat();
      }
      tswitch = R::rbinom(1, p(da_cat));
    }
  else if (model == "das28-switch"){ // Treatment -> DA -> Switch
        das28 = das28 + nma_dhaq.sim_dy();
        da_cat = get_das28_cat();
        tswitch = R::rbinom(1, p(da_cat));
  }
  else if (model == "acr-eular-switch"){ // Treatment -> ACR -> EULAR -> Switch
    if (eular == 0){
      tswitch = 1;
    }
  }
  
  // NBT
  if (therapy == nbt){
    tswitch = 0;
  }
}

/***********************************
* Time to treatment discontinuation
***********************************/
// Time to treatment discontinuation by eular response 
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

/**********************************
* Latent class growth model (LCGM)
**********************************/
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

/****************
* Simulate costs
****************/
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

//// Hospitalization costs and days
struct Hosp {
  double cost;
  double days;
};

Hosp sim_hosp(double &haq, double yrlen, arma::rowvec hosp_days, 
                       arma::rowvec cost_pday){
  Hosp hosp;
  if (haq < 0.5){
    hosp.days = hosp_days(0) * yrlen;
    hosp.cost = hosp.days * cost_pday(0);
  }
  else if (haq >= 0.5 && haq < 1){
    hosp.days = hosp_days(1) * yrlen;
    hosp.cost = hosp.days * cost_pday(1);
  }
  else if (haq >= 1 && haq < 1.5){
    hosp.days = hosp_days(2) * yrlen;
    hosp.cost = hosp.days * cost_pday(2);
  }
  else if (haq >= 1.5 & haq < 2){
    hosp.days = hosp_days(3) * yrlen;
    hosp.cost = hosp.days * cost_pday(3);
  } 
  else if (haq >= 2 & haq < 2.5){
    hosp.days = hosp_days(4) * yrlen;
    hosp.cost = hosp.days * cost_pday(4);
  }
  else if (haq >= 2.5){
    hosp.days = hosp_days(5) * yrlen;
    hosp.cost = hosp.days * cost_pday(5);
  }
  return hosp;
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

/****************************
* Simulate utility and QALYs
****************************/
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

/**************************************************************
* Calculate summary statistics during simulation run
* Note that this is particularly usefeul when the user wants 
* to save RAM and would rather have summary results than
* output for each individual and model cycle
**************************************************************/
// Calculate means by individual for selected variables
class SimMeans {         
private:
  int mod;
  int sim; 
  int indiv;
  int cycle;
  int index;
  int n_indivs;
  int n_sims;
  int n_mods;
  double cycle_length;
  double yrlen;
  double discount_qalys;
  double discount_cost;
  double beta_qalys;
  double beta_cost;
  std::map<std::string, double > indivsums;
  std::map<std::string, std::vector<int> > id;
  std::map<std::string, std::vector<double> > varsums;
public:  
  SimMeans(int n_mods, int n_sims_, int n_indivs_, double r_qalys = .03, double r_cost = .03,
           double cycle_legth_ = 6);
  std::map<std::string, std::vector<double> > get_varsums();
  std::map<std::string, std::vector<int> > get_id();
  void set_iterators(int m, int s, int i, int c);
  void set_id();
  void increment_id(int m, int s, int i, int c);
  double calc_lys_infusion(std::string route, double &lys);
  double calc_lys_injection(std::string route, double &lys);
  double calc_lys_oral(std::string route, double &lys);
  double calc_dhaq(double haq0, double haq, bool final_cycle);
  void increment_indivsums(double yrs_since_approval, double dqalys, double dhc_cost, 
                           double dprod_loss);
  void increment_varsums(double qalys, double tx_cost, Hosp hosp,
                         double mgmt_cost, double si_cost, double prod_loss, double si,
                         std::string route, double haq0, double haq, bool final_cycle,
                         double yrs_since_approval, double age0); 
  std::map<std::string, std::vector<double> > calc_means();
};

SimMeans::SimMeans(int n_mods_, int n_sims_, int n_indivs_, double r_qalys, double r_cost,
                   double cycle_length_){
  cycle = 0;
  indiv = 0;
  sim = 0;
  mod = 0;
  cycle_length = cycle_length_;
  yrlen = cycle_length/12;
  n_indivs = n_indivs_;
  n_sims = n_sims_;
  n_mods = n_mods_;
  index = 0;
  discount_qalys = r_qalys;
  discount_cost = r_cost;
  beta_qalys = pow(1/(1 + discount_qalys), cycle_length/12); 
  beta_cost = pow(1/(1 + discount_cost), cycle_length/12); 
  int N = n_sims * n_mods;
  indivsums = std::map<std::string,double >();
  id = std::map<std::string, std::vector<int> > ();
  varsums = std::map<std::string, std::vector<double> >();
  indivsums["yrs_since_approval"] = 0.0;
  indivsums["dqalys"] = 0.0;
  indivsums["dhc_cost"] = 0.0;
  indivsums["dprod_loss"] = 0.0;
  id["sim"] = std::vector<int> (N);
  id["mod"] = std::vector<int> (N);
  varsums["lys"] = std::vector<double> (N);
  varsums["dlys"] = std::vector<double> (N);
  varsums["lys_infusion"] = std::vector<double> (N);
  varsums["lys_injection"] = std::vector<double> (N);
  varsums["lys_oral"] = std::vector<double> (N);
  varsums["dhaq"] = std::vector<double> (N);
  varsums["si"] = std::vector<double> (N);
  varsums["qalys"] = std::vector<double> (N);
  varsums["dqalys"] = std::vector<double> (N);
  varsums["tx_cost"] = std::vector<double> (N);
  varsums["dtx_cost"] = std::vector<double> (N);
  varsums["hosp_days"] = std::vector<double> (N);
  varsums["hosp_cost"] = std::vector<double> (N);
  varsums["dhosp_cost"] = std::vector<double> (N);
  varsums["mgmt_cost"] = std::vector<double> (N);
  varsums["dmgmt_cost"] = std::vector<double> (N);
  varsums["si_cost"] = std::vector<double> (N);
  varsums["dsi_cost"] = std::vector<double> (N);
  varsums["prod_loss"] = std::vector<double> (N);
  varsums["dprod_loss"] = std::vector<double> (N);
  varsums["dhc_cost"] = std::vector<double> (N);
  varsums["dtot_cost"] = std::vector<double> (N);
  varsums["yrs_since_approval"] = std::vector<double> (N);
  varsums["dqalys_ann"] = std::vector<double> (N);
  varsums["dhc_cost_ann"] = std::vector<double> (N);
  varsums["dprod_loss_ann"] = std::vector<double> (N);
}

std::map<std::string, std::vector<double> > SimMeans::get_varsums(){
  return varsums;
}

std::map<std::string, std::vector<int> > SimMeans::get_id(){
  return id;
}

void SimMeans::set_iterators(int m, int s, int i, int c){
  mod = m;
  sim = s;
  indiv = i;
  cycle = c;
  index = mod * n_sims + sim;
}

void SimMeans::set_id(){
  id["mod"][index] = mod;
  id["sim"][index] = sim;
}

void SimMeans::increment_id(int m, int s, int i, int c){
  set_iterators(m, s, i, c);
  set_id();
}

double  SimMeans::calc_lys_infusion(std::string route, double &lys){
  if (route == "infusion"){
    return lys;
  }
  else if (route == "infusion/injection"){
    return lys/2;
  } 
  else{
    return 0;
  }
}

double SimMeans::calc_lys_injection(std::string route, double &lys){
  if (route == "injection"){
    return lys;
  }
  else if (route == "infusion/injection" || route == "oral/injection"){
    return lys/2;
  } 
  else{
    return 0;
  }
}

double SimMeans::calc_lys_oral(std::string route, double &lys){
  if (route == "oral"){
    return lys;
  }
  else if (route == "oral/injection"){
    return lys/2;
  } 
  else{
    return 0;
  }
}

double SimMeans::calc_dhaq(double haq0, double haq, bool final_cycle){
  if (final_cycle == false){
    return 0;
  }
  else{
    return haq - haq0;
  }
}

void SimMeans::increment_indivsums(double yrs_since_approval, double dqalys, double dhc_cost,
                                   double dprod_loss){
  if (cycle == 0){
    indivsums["yrs_since_approval"] = 0.0 + yrs_since_approval;
    indivsums["dqalys"] = 0.0 + dqalys;
    indivsums["dhc_cost"] = 0.0 + dhc_cost;
    indivsums["dprod_loss"] = 0.0 + dprod_loss;
  }
  else{
    indivsums["yrs_since_approval"] = indivsums["yrs_since_approval"] + yrs_since_approval;
    indivsums["dqalys"] = indivsums["dqalys"] + dqalys;
    indivsums["dhc_cost"] = indivsums["dhc_cost"] + dhc_cost;
    indivsums["dprod_loss"] = indivsums["dprod_loss"] + dprod_loss;
  }
  
}

void SimMeans::increment_varsums(double qalys, double tx_cost, Hosp hosp,
                                 double mgmt_cost, double si_cost, double prod_loss, double si,
                                 std::string route, double haq0, double haq, bool final_cycle,
                                 double yrs_since_approval, double age0){
  double dfq = discount_factor(cycle + 1, discount_qalys, yrlen);
  double dfc = discount_factor(cycle + 1, discount_cost, yrlen);
  double dqalys = qalys * dfq;
  double dtx_cost = tx_cost * dfc;
  double dhosp_cost = hosp.cost * dfc;
  double dmgmt_cost = mgmt_cost * dfc;
  double dsi_cost = si_cost * dfc;
  double dprod_loss = prod_loss * dfc;
  double dhc_cost = dtx_cost + dhosp_cost + dmgmt_cost + dsi_cost;
  double dtot_cost = dhc_cost + dprod_loss;
  increment_indivsums(yrs_since_approval, dqalys, dhc_cost, dprod_loss);
  varsums["lys"][index] = varsums["lys"][index] + yrlen;
  varsums["dlys"][index] = varsums["dlys"][index] + yrlen * dfq;
  varsums["lys_infusion"][index] = varsums["lys_infusion"][index] + calc_lys_infusion(route, yrlen);
  varsums["lys_injection"][index] = varsums["lys_injection"][index] + calc_lys_injection(route, yrlen);
  varsums["lys_oral"][index] = varsums["lys_oral"][index] + calc_lys_oral(route, yrlen);
  varsums["dhaq"][index] = varsums["dhaq"][index] + calc_dhaq(haq0, haq, final_cycle);
  varsums["si"][index] = varsums["si"][index] + si;
  varsums["qalys"][index] = varsums["qalys"][index] + qalys;
  varsums["dqalys"][index] = varsums["dqalys"][index] + dqalys;
  varsums["tx_cost"][index] = varsums["tx_cost"][index] + tx_cost;
  varsums["dtx_cost"][index] = varsums["dtx_cost"][index] + dtx_cost;
  varsums["hosp_days"][index] = varsums["hosp_days"][index] + hosp.days;
  varsums["hosp_cost"][index] = varsums["hosp_cost"][index] + hosp.cost;
  varsums["dhosp_cost"][index] = varsums["dhosp_cost"][index] + dhosp_cost;
  varsums["mgmt_cost"][index] = varsums["mgmt_cost"][index] + mgmt_cost;
  varsums["dmgmt_cost"][index] = varsums["dmgmt_cost"][index] + dmgmt_cost;
  varsums["si_cost"][index] = varsums["si_cost"][index] + si_cost;
  varsums["dsi_cost"][index] = varsums["dsi_cost"][index] + dsi_cost;
  varsums["prod_loss"][index] = varsums["prod_loss"][index] + prod_loss;
  varsums["dprod_loss"][index] = varsums["dprod_loss"][index] + dprod_loss;
  varsums["dhc_cost"][index] = varsums["dhc_cost"][index] + dhc_cost;
  varsums["dtot_cost"][index] = varsums["dtot_cost"][index] + dtot_cost;
  if (final_cycle == true){
    //int maxt = (100 - int(age0)) * (12/cycle_length);
    varsums["yrs_since_approval"][index] = varsums["yrs_since_approval"][index] +
      indivsums["yrs_since_approval"]/(cycle + 1);
    varsums["dqalys_ann"][index] = varsums["dqalys_ann"][index] +
      const_geometric_series(indivsums["dqalys"], beta_qalys, cycle + 1, 1) * 12/cycle_length;
    varsums["dhc_cost_ann"][index] = varsums["dhc_cost_ann"][index] +
      const_geometric_series(indivsums["dhc_cost"], beta_cost, cycle + 1, 1) * 12/cycle_length;
    varsums["dprod_loss_ann"][index] = varsums["dprod_loss_ann"][index] +
      const_geometric_series(indivsums["dprod_loss"], beta_cost, cycle + 1, 1) * 12/cycle_length;
  }
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
  void increment_varsums(double qalys, double haq, double tx_cost, double hosp_cost, double mgmt_cost, 
                         double si_cost, double prod_loss); 
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
  varsums["tx_cost"] = std::vector<double> (N);
  varsums["hosp_cost"] = std::vector<double> (N);
  varsums["mgmt_cost"] = std::vector<double> (N);
  varsums["si_cost"] = std::vector<double> (N);
  varsums["prod_loss"] = std::vector<double> (N);
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

void TimeMeans::increment_varsums(double qalys, double haq, double tx_cost, double hosp_cost, double mgmt_cost, 
                                  double si_cost, double prod_loss){
  varsums["qalys"][index] = varsums["qalys"][index] + qalys;
  varsums["haq"][index] = varsums["haq"][index] + haq;
  varsums["tx_cost"][index] = varsums["tx_cost"][index] + tx_cost;
  varsums["hosp_cost"][index] = varsums["hosp_cost"][index] + hosp_cost;
  varsums["mgmt_cost"][index] = varsums["mgmt_cost"][index] + mgmt_cost;
  varsums["si_cost"][index] = varsums["si_cost"][index] + si_cost;
  varsums["prod_loss"][index] = varsums["prod_loss"][index] + prod_loss;
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
// [[Rcpp::export]]
std::vector<double> sim_qalysC(std::vector<double> &utility, std::vector<double> &yrlen,
                           std::vector<int> &sim, std::vector<int> &tx,
                           std::vector<int> &si, std::vector<double> &si_ul,
                           arma::mat &x_attr, arma::mat &tx_attr_coef){
  int N = utility.size();
  std::vector<double> qalys_vec;
  qalys_vec.reserve(N);
  for (int i = 0; i < N; ++i){
    double utility1 = std::min(1.0, utility[i] - si[i] * si_ul[sim[i]]/12 + 
      arma::dot(x_attr.row(tx[i]), tx_attr_coef.row(sim[i])));
    qalys_vec.push_back(yrlen[i] * utility1);
  }
  return qalys_vec;
}

/*********************************
* The MAIN C++ sim_iviRA function
*********************************/
// Simulate HAQ score
// [[Rcpp::export]] 
List sim_iviRA_C(arma::mat arm_inds, Rcpp::DataFrame tx_data, 
                 CharacterMatrix model_structures_mat, std::string hist,
             std::vector<double> haq0, std::vector<double> das28_0,
             std::vector<double> sdai0, std::vector<double> cdai0,
             std::vector<double> age0, std::vector<int> male,
             std::vector<int> prev_dmards,
             Rcpp::List nma_acr_list, arma::mat x_acr,
             Rcpp::List nma_haq_list, arma::mat x_haq,
             Rcpp::List nma_das28_list, arma::mat x_das28,
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
             Rcpp::List utility_tx_attr,
             Rcpp::List discount_rate, std::string output){
  
  // Type conversions
  //// Treatment data
  std::vector<std::string> route = as<std::vector<std::string> >(tx_data["route"]);
  std::vector<double> yrs_since_approval = as<std::vector<double> >(tx_data["years_since_approval"]);
  
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
  
  // NMA ACR response
  std::vector<double> acr_rr = as<std::vector<double> > (nma_acr_list["rr"]);
  std::vector<double> acr_A = as<std::vector<double> > (nma_acr_list["A"]);
  std::vector<double> acr_z2 = as<std::vector<double> > (nma_acr_list["z2"]);
  std::vector<double> acr_z3 = as<std::vector<double> > (nma_acr_list["z3"]);
  arma::cube acr_d_beta = as<arma::cube>(nma_acr_list["d"]);
  
  // NMA HAQ
  std::vector<double> nma_haq_rr = as<std::vector<double> > (nma_haq_list["rr"]);
  std::vector<double> nma_haq_A = as<std::vector<double> > (nma_haq_list["A"]);
  arma::cube nma_haq_d_beta = as<arma::cube>(nma_haq_list["d"]);
  
  // NMA DAS28
  std::vector<double> nma_das28_rr = as<std::vector<double> > (nma_das28_list["rr"]);
  std::vector<double> nma_das28_A = as<std::vector<double> > (nma_das28_list["A"]);
  arma::cube nma_das28_d_beta = as<arma::cube>(nma_das28_list["d"]);
  
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
  
  //// Treatment attributes
  arma::mat x_attr = as<arma::mat> (utility_tx_attr["data"]);
  arma::mat tx_attr_utilpars = as<arma::mat> (utility_tx_attr["pars"]);
  
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
  nmaACR nma_acr;
  nmaLM nma_haq;
  nmaLM nma_das28;
  TxIHaq tx_ihaq;
  Cost sim_cost;
  ModelStructure mod_struct;
  double utility = 0.0;
  arma::rowvec utilmix_x(5);
  arma::rowvec utilmix_w(4);
  double qalys = 0.0;
  SimMeans sim_means(n_mods, n_sims, n_pat, discount_qalys, discount_cost, cycle_length);
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
  std::vector<double> hosp_days_vec;
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
        double month = 6; // This is because RCTs are 6 months long!
        cage_i = age0[i];
        age_i = int(cage_i);
        double haq = haq0[i];
        int cycle = 0;
        bool final_cycle = false;
        if (arm_inds.n_rows > 1){
          arm_row = i;
          arm_ind_i = arm_inds.row(i);
        }
        int t_cdmards = 0;
        TxISwitch tx_iswitch;
        tx_iswitch.das28 = das28_0[i];
        tx_iswitch.sdai = sdai0[i];
        tx_iswitch.cdai = cdai0[i];
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
          arma::rowvec x_attr_ij = x_attr.row(arm_ind_ij);
          
          // H1-H3: simulate change in HAQ during initial treatment phase
          nma_acr.set(hist, acr_rr[s], acr_A[s], acr_z2[s], acr_z3[s],
                      acr_d_beta.slice(arm_ind_ij).row(s), x_acr.row(i),
                      j);
          nma_haq.set(hist, nma_haq_rr[s], nma_haq_A[s],
                      nma_haq_d_beta.slice(arm_ind_ij).row(s), x_haq.row(i),
                      j);
          nma_das28.set(hist, nma_das28_rr[s], nma_das28_A[s],
                      nma_das28_d_beta.slice(arm_ind_ij).row(s), x_das28.row(i),
                      j);
          tx_ihaq.sim(mod_struct.tx_ihaq, j, arm_ind_ij, nbt,
                                          nma_acr, nma_haq,
                                          acr2eular.slice(s), acr2haq.row(s), 
                                          eular2haq.row(s));
          
          // S1-S6: simulate treatment switching during first 6 months
          tx_iswitch.sim(mod_struct.tx_iswitch, j, arm_ind_ij, nbt,
                        tx_ihaq.acr, tx_ihaq.eular,
                        acr2das28.row(s), acr2sdai.row(s), acr2cdai.row(s),
                        nma_das28, tswitch_da.row(s));
          
          // Time to treatment discontinuation and serious infections
          double ttsi_j = (-cycle_length/12 + rsurvC(as_scalar(si_loc.row(s).col(arm_ind_ij)), 
                                                as_scalar(si_anc1.row(s).col(arm_ind_ij)),
                                                si_dist,as_scalar(si_anc2.row(s).col(arm_ind_ij))
                                                )) * (12/cycle_length);
          double ttd_j = 0;
          if (mod_struct.ttd_cause == "all"){
            if (mod_struct.tx_iswitch == "acr-eular-switch"){
              ttd_j = sim_ttd_eular(x_ttd_eular.row(i), 
                                    ttd_eular_mod.loc[mod_struct.ttd_dist].row(s), ttd_eular_mod.anc1[mod_struct.ttd_dist](s), 
                                    ttd_eular_good.loc[mod_struct.ttd_dist].row(s), ttd_eular_good.anc1[mod_struct.ttd_dist](s),
                                    tx_ihaq.eular, mod_struct.ttd_dist, cycle_length, ttsi_j,
                                    ttd_eular_mod.anc2[mod_struct.ttd_dist](s), ttd_eular_good.anc2[mod_struct.ttd_dist](s));
            }
            else if (mod_struct.tx_iswitch == "acr-switch"){
              ttd_j = sim_ttd(x_ttd_all.row(i), ttd_all.loc[mod_struct.ttd_dist].row(s), ttd_all.anc1[mod_struct.ttd_dist](s),
                              tx_iswitch.tswitch, mod_struct.ttd_dist, cycle_length, ttsi_j,
                              ttd_all.anc2[mod_struct.ttd_dist](s));
            }
            else {
              x_ttd_da_i = update_x_ttd_da(x_ttd_da.row(i), tx_iswitch.da_cat);
              ttd_j = sim_ttd(x_ttd_da.row(i), ttd_da.loc[mod_struct.ttd_dist].row(s), ttd_da.anc1[mod_struct.ttd_dist](s),
                              tx_iswitch.tswitch, mod_struct.ttd_dist, cycle_length, ttsi_j,
                              ttd_da.anc2[mod_struct.ttd_dist](s));
            }  
          }
          else if (mod_struct.ttd_cause == "si"){
            if (ttsi_j < 0){
                ttd_j = 0;
            } else {
                ttd_j = ttsi_j;
            }
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
            int death = sample_deathC(age_i, male[i], lifetable_male, lifetable_female, 
                                         x_mort.row(i), logor_mort.row(s), haq0[i], haq,
                                        cycle_length, month, haqdelta_loghr.row(s));
            if (death == 1 || max_months < month + cycle_length){
              final_cycle = true;
            }
            
            // Did serious infection occur?
            if ((ttsi_j < ttd_j && t > ttd_j && t < ttd_j + 1 && arm_ind_ij != nbt) ||
                (ttsi_j < 0 && arm_ind_ij != nbt) ||
                (mod_struct.ttd_cause == "si" && t > ttd_j && t < ttd_j + 1 &&
                arm_ind_ij != nbt)) {
                  si = 1;
            } 
            
            // Update HAQ score
            if (t == 0){ // initial haq change
              if (ttd_j == 0){ // immediate rebound for patients switching treatment during the initial period
                  // haq = haq + tx_ihaq.dhaq - tx_ihaq.dhaq * rebound_factor[s]; 
                  haq = haq;
              } else{
                haq = haq + tx_ihaq.dhaq; 
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
              haq = haq - tx_ihaq.dhaq * rebound_factor[s];
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
            Hosp hosp = sim_hosp(haq, cycle_length/12, hosp_days.row(s), cost_pday.row(s));
            sim_cost.hosp = hosp.cost;
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
            utility = std::min(1.0, utility - si * si_ul[s]/12 + 
              arma::dot(x_attr_ij, tx_attr_utilpars.row(s)));
            qalys = cycle_length/12 * utility;
            
            if (output == "summary"){
              sim_means.increment_id(m, s, i, cycle);
              sim_means.increment_varsums(qalys, sim_cost.tx, hosp,
                                          sim_cost.mgmt, sim_cost.si, sim_cost.prod, si,
                                          route[arm_ind_ij], haq0[i], haq, final_cycle, 
                                          yrs_since_approval[arm_ind_ij], age0[i]);
              time_means.increment_id(m, s, month); //note: requires assumption of constant model cycles!!!
              time_means.increment_alive();
              time_means.increment_varsums(qalys, haq, sim_cost.tx, sim_cost.hosp, sim_cost.mgmt, 
                                           sim_cost.si, sim_cost.prod);
              out0.push_back(t, m, s, i, arm_ind_ij, tx_ihaq.acr, tx_ihaq.eular,
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
            
            acr_vec.push_back(tx_ihaq.acr);
            eular_vec.push_back(tx_ihaq.eular);
            das28_vec.push_back(tx_iswitch.das28);
            sdai_vec.push_back(tx_iswitch.sdai);
            cdai_vec.push_back(tx_iswitch.cdai);
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
            hosp_days_vec.push_back(hosp.days);
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
            ++cycle;
            
            // Stop i loop for death or final model cycle
            if (final_cycle == true){
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
      _["hosp_days"] = hosp_days_vec,
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
    Rcpp::DataFrame sim_means_df1 = Rcpp::DataFrame::create(
      _["model"] = sim_means_id["mod"],
      _["sim"] = sim_means_id["sim"],
      _["lys"] = sim_means_out["lys"],
      _["dlys"] = sim_means_out["dlys"],
      _["lys_infusion"] = sim_means_out["lys_infusion"],
      _["lys_injection"] = sim_means_out["lys_injection"],
      _["lys_oral"] = sim_means_out["lys_oral"],
      _["dhaq"] = sim_means_out["dhaq"],                             
      _["si"] = sim_means_out["si"],                        
      _["qalys"] = sim_means_out["qalys"],
      _["dqalys"] = sim_means_out["dqalys"], 
      _["tx_cost"] = sim_means_out["tx_cost"],
      _["dtx_cost"] = sim_means_out["dtx_cost"], 
      _["hosp_days"] = sim_means_out["hosp_days"],
      _["hosp_cost"] = sim_means_out["hosp_cost"],                             
      _["dhosp_cost"] = sim_means_out["dhosp_cost"],
      _["mgmt_cost"] = sim_means_out["mgmt_cost"],
      _["dmgmt_cost"] = sim_means_out["dmgmt_cost"],
      _["si_cost"] = sim_means_out["si_cost"],
      _["dsi_cost"] = sim_means_out["dsi_cost"]
    );
    Rcpp::DataFrame sim_means_df2 = Rcpp::DataFrame::create(
      _["prod_loss"] = sim_means_out["prod_loss"],
      _["dprod_loss"] = sim_means_out["dprod_loss"],
      _["dhc_cost"] = sim_means_out["dhc_cost"],
      _["dtot_cost"] = sim_means_out["dtot_cost"],
      _["yrs_since_approval"] = sim_means_out["yrs_since_approval"],
      _["dqalys_ann"] = sim_means_out["dqalys_ann"],
      _["dhc_cost_ann"] = sim_means_out["dhc_cost_ann"],
      _["dprod_loss_ann"] = sim_means_out["dprod_loss_ann"]
    );
    std::map<std::string, std::vector<double> > time_means_out = time_means.calc_means();
    std::map<std::string, std::vector<int> > time_means_id = time_means.get_id();
    Rcpp::DataFrame time_means_df = Rcpp::DataFrame::create(
      _["model"] = time_means_id["mod"],
      _["sim"] = time_means_id["sim"],
      _["month"] = time_means_id["month"],
      _["alive"] = time_means.get_alive(),
      _["qalys"] = time_means_out["qalys"],
      _["haq"] = time_means_out["haq"],
      _["tx_cost"] = time_means_out["tx_cost"],
      _["hosp_cost"] = time_means_out["hosp_cost"],
      _["mgmt_cost"] = time_means_out["mgmt_cost"],
      _["si_cost"] = time_means_out["si_cost"],
      _["prod_loss"] = time_means_out["prod_loss"]
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
      _["means1"] = sim_means_df1,
      _["means2"] = sim_means_df2,
      _["time.means"] = time_means_df,
      _["out0"] = out0_df
      );
  }
}



