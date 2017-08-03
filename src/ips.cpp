// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(hesim)]]

#include <RcppArmadillo.h>
#include "rand.h"
#include "mort.h"
#include "logit-prob.h"
#include <hesim.h>
using namespace Rcpp;

// Update HAQ initial response
void update_haq_t1(double &haq, double haq_change){
  haq = haq + haq_change;
}

// Effect of treatment on HAQ during the initial treatment phase
struct ItreatHaq {
  int acr;
  int eular;
  double dhaq;
};

ItreatHaq sim_itreat_haq(std::string itreat_haq_model, int line, int therapy, int nbt,
                     arma::rowvec nma_acr1, arma::rowvec nma_acr2, 
                     double nma_dhaq1, double nma_dhaq2,
                     arma::mat acr2eular, arma::rowvec acr2haq, arma::rowvec eular2haq){
  ItreatHaq sim;
  
  // ACR response
  if (itreat_haq_model == "acr-haq" || itreat_haq_model == "acr-eular-haq"){
      if (line == 0){
          sim.acr = hesim::rcat1C(nma_acr1);
        }  
      else{
          sim.acr = hesim::rcat1C(nma_acr2);
      }
  }
  
  // EULAR response
  if (itreat_haq_model == "acr-eular-haq"){
      sim.eular = hesim::rcat1C(acr2eular.row(sim.acr));
  }
  
  // HAQ
  if (itreat_haq_model == "acr-haq"){
      sim.dhaq = acr2haq(sim.acr);
  } 
  else if (itreat_haq_model == "acr-eular-haq"){
      sim.dhaq = eular2haq(sim.eular);
  } 
  else if (itreat_haq_model == "haq") {
      if (line == 0){
        sim.dhaq = nma_dhaq1;
      } else {
        sim.dhaq = nma_dhaq2;
      }
  }
  
  // No treatment response for nbt
  if (therapy == nbt){
    sim.acr = 0;
    sim.eular = 0;
    sim.dhaq = 0;
  }
  return sim;
}

// [[Rcpp::export]]
List test_itreat_haq(std::string itreat_haq_model, int line, int therapy, int nbt,
                      arma::rowvec nma_acr1, arma::rowvec nma_acr2, 
                      double nma_dhaq1, double nma_dhaq2,
                      arma::mat acr2eular, arma::rowvec acr2haq, arma::rowvec eular2haq){
  ItreatHaq sim = sim_itreat_haq(itreat_haq_model, line, therapy, nbt,
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
struct ItreatSwitch {
  int tswitch;
  int das28;
  int sdai;
  int cdai;
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

ItreatSwitch sim_itreat_switch(std::string itreat_switch_model, int line, int therapy, int nbt,
                  int acr, int eular, double das28, double sdai, double cdai,
                  arma::rowvec acr2das28, arma::rowvec acr2sdai,arma::rowvec acr2cdai,
                  arma::rowvec nma_das28_1, arma::rowvec nma_das28_2, arma::rowvec p){
  
  ItreatSwitch sim;
  sim.tswitch = 0;
  sim.da_cat = 0;
  int da_change = 0;
  
  // Switching rules
  if (itreat_switch_model == "acr-switch"){ // Treatment -> ACR -> Switch
    if (acr == 0){
      sim.tswitch = 1;
    }
  }
  else if (itreat_switch_model == "acr-das28-switch" ||   // Treatment -> ACR -> DA -> Switch
           itreat_switch_model == "acr-sdai-switch" ||
           itreat_switch_model == "acr-cdai-switch") {
      if (itreat_switch_model == "acr-das28-switch"){
        da_change = acr2das28(acr);
        double das28_new = get_da_new(das28, da_change, 0, 9.4);
        sim.da_cat = get_das28_cat(das28_new);
      }
      else if (itreat_switch_model == "acr-sdai-switch"){
        da_change = acr2sdai(acr);
        double sdai_new = get_da_new(sdai, da_change, 0, 86);
        sim.da_cat = get_sdai_cat(sdai_new);
      }
      else if (itreat_switch_model == "acr-cdai-switch"){
        da_change = acr2cdai(acr);
        double cdai_new = get_da_new(cdai, da_change, 0, 76);
        sim.da_cat = get_cdai_cat(cdai_new);
      }
      sim.tswitch = R::rbinom(1, p(sim.da_cat));
    }
  else if (itreat_switch_model == "das28-switch"){ // Treatment -> DA -> Switch
      if (line == 0){
        sim.da_cat = get_das28_cat(das28 + nma_das28_1(therapy));
      }
      else{
        sim.da_cat = get_das28_cat(das28 + nma_das28_2(therapy));
      }
      sim.tswitch = R::rbinom(1, p(sim.da_cat));
  }
  else if (itreat_switch_model == "acr-eular-switch"){ // Treatment -> ACR -> EULAR -> Switch
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
List test_itreat_switch(std::string itreat_switch_model, int line, int therapy, int nbt,
                        int acr, int eular, double das28, double sdai, double cdai,
                        arma::rowvec acr2das28, arma::rowvec acr2sdai, arma::rowvec acr2cdai,
                        arma::rowvec nma_das28_1, arma::rowvec nma_das28_2, arma::rowvec p){
  ItreatSwitch sim = sim_itreat_switch(itreat_switch_model, line, therapy, nbt, acr, eular,
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
                     double discount){
  
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
          dose_num_j * price[agents_ind_j] * (1 - discount) + infusion_cost_j;
      }
      else{
        tc += (ceil(weight * dose_amt[agents_ind_j]/strength_val[agents_ind_j]) * 
          dose_num_j * price[agents_ind_j] * (1 - discount) + 
          infusion_cost[agents_ind_j] * dose_num[agents_ind_j]) * cycle_length/12; 
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

// Simulate HAQ score
// [[Rcpp::export]] 
List sim_iviRA_C(arma::mat arm_inds,
             std::vector<double> haq0, std::vector<double> das28_0,
             std::vector<double> sdai0, std::vector<double> cdai0,
             std::vector<double> age0, std::vector<int> male,
             std::vector<int> prev_dmards,
             std::string itreat_haq_model, std::string itreat_switch_model,
             arma::cube nma_acr1, arma::cube nma_acr2, arma::mat nma_dhaq1, arma::mat nma_dhaq2,
             arma::mat nma_das28_1, arma::mat nma_das28_2,
             arma::cube acr2eular,arma::mat acr2haq, arma::mat eular2haq, 
             arma::mat acr2das28, arma::mat acr2sdai, arma::mat acr2cdai,
             arma::mat tswitch_da,
             arma::mat haq_lprog_therapy, arma::mat haq_lprog_age,
             arma::cube haq_lcgm_delta, arma::cube haq_lcgm_beta, std::string cdmards_haq_model,
             std::vector<double> rebound_factor,
             arma::mat lifetable_male, arma::mat lifetable_female, 
             arma::mat x_mort, arma::mat logor_mort, 
             std::string ttd_dist, arma::mat x_ttd, 
             arma::mat ttd_all_loc, arma::vec ttd_all_anc1, arma::vec ttd_all_anc2,
             arma::mat ttd_da_loc, arma::vec ttd_da_anc1, arma::vec ttd_da_anc2,
             arma::mat ttd_eular_loc_mod, arma::vec ttd_eular_anc1_mod, arma::vec ttd_eular_anc2_mod,
             arma::mat ttd_eular_loc_good, arma::vec ttd_eular_anc1_good, arma::vec ttd_eular_anc2_good,
             double cycle_length, int treat_gap, int cdmards, int nbt, 
             arma::mat si_loc, arma::mat si_anc1, arma::mat si_anc2, std::string si_dist, 
             arma::mat haqdelta_loghr, int max_months, 
             arma::mat hosp_days, arma::mat cost_pday, std::vector<double> mgmt_cost, 
             std::vector<double> si_cost, std::vector<double> prod_loss, 
             arma::cube tc_agents_ind, Rcpp::List tc_list, std::vector<double> weight){
  
  // Type conversions
  std::vector<double> tc_discount = as<std::vector<double> >(tc_list["discount"]);
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

  // Constants
  int n_treatments = arm_inds.n_cols;
  int n_sims = logor_mort.n_rows;
  int n_pat = x_mort.n_rows;
  int maxt = 100 * (12/cycle_length);
  
  // Declarations
  int counter = 0;
  double cage_i = 0.0; // continuous age
  int age_i = 0;
  int arm_row = 0;
  arma::rowvec arm_ind_i = arm_inds.row(0);
  int arm_ind_ij = 0;
  Cost sim_cost;
  
  // Information to store
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
  std::vector<double> haq_vec;
  std::vector<double> ttsi_vec;
  std::vector<int> si_vec;
  std::vector<double> yrlen_vec;
  std::vector<double> tx_cost_vec;
  std::vector<double> hosp_cost_vec;
  std::vector<double> mgmt_cost_vec;
  std::vector<double> si_cost_vec;
  std::vector<double> prod_loss_vec;
  
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
      ItreatSwitch sim_s_t1;
      sim_s_t1.das28 = das28_0[i];
      sim_s_t1.sdai = sdai0[i];
      sim_s_t1.cdai = cdai0[i];
      
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
        ItreatHaq sim_h_t1 = sim_itreat_haq(itreat_haq_model, j, arm_ind_ij, nbt,
                                        nma_acr1.slice(arm_ind_ij).row(s), 
                                        nma_acr2.slice(arm_ind_ij).row(s), 
                                        nma_dhaq1.col(arm_ind_ij)(s),
                                        nma_dhaq2.col(arm_ind_ij)(s),
                                        acr2eular.slice(s), acr2haq.row(s), 
                                        eular2haq.row(s));
        
        // S1-S6: simulate treatment switching during first 6 months
        sim_s_t1 = sim_itreat_switch(itreat_switch_model, j, arm_ind_ij, nbt,
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
        if (itreat_switch_model == "acr-eular-switch"){
            ttd_j = sim_ttd_eular(x_ttd.row(i), ttd_eular_loc_mod.row(s), ttd_eular_anc1_mod(s), 
                                                 ttd_eular_loc_good.row(s), ttd_eular_anc1_good(s), 
                                                 sim_h_t1.eular, ttd_dist, cycle_length, ttsi_j,
                                                 ttd_eular_anc2_mod(s), ttd_eular_anc2_good(s));
        }
        else if (itreat_switch_model == "acr-switch"){
          ttd_j = sim_ttd(x_ttd.row(i), ttd_all_loc.row(s), ttd_all_anc1(s),
                              sim_s_t1.tswitch, ttd_dist, cycle_length, ttsi_j,
                              ttd_all_anc2(s));
        }
        else {
          ttd_j = sim_ttd(x_ttd.row(i), ttd_da_loc.row(s), ttd_da_anc1(s),
                          sim_s_t1.tswitch, ttd_dist, cycle_length, ttsi_j,
                          ttd_all_anc2(s));
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
            if(cdmards_haq_model == "lcgm" && (arm_ind_ij == nbt || arm_ind_ij == cdmards)){
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
          
          // Update costs
          sim_cost.tx = sim_tx_cost1C(t, tc_agents_ind.slice(arm_row).row(j), tx_name,
                                      init_dose_val, ann_dose_val, strength_val,
                                      init_num_doses, ann_num_doses, price,
                                      infusion_cost, loading_dose, 
                                      weight_based, weight[i], cycle_length, tc_discount[s]);
          sim_cost.hosp = sim_hosp_cost1C(haq, cycle_length/12, hosp_days.row(s), cost_pday.row(s));
          sim_cost.mgmt = sim_mgmt_cost1C(cycle_length/12, mgmt_cost[s]);
          sim_cost.prod = sim_prod_loss1C(haq, cycle_length/12, prod_loss[s]);
          sim_cost.si = sim_si_cost1C(si, cycle_length/12, si_cost[s]);
          
          // Fill vectors
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
  return List::create(Rcpp::Named("sim") = sim, Rcpp::Named("id") = id, 
                      Rcpp::Named("month") = month_vec, Rcpp::Named("tx") = tx_vec,
                      Rcpp::Named("tx_seq") = tx_seq, Rcpp::Named("tx_cycle") = tx_cycle,
                      Rcpp::Named("age") = age_vec, Rcpp::Named("ttd") = ttd_vec, 
                      Rcpp::Named("acr") = acr_vec, Rcpp::Named("eular") = eular_vec,
                      Rcpp::Named("haq") = haq_vec, Rcpp::Named("ttsi") = ttsi_vec, 
                      Rcpp::Named("si") = si_vec, Rcpp::Named("yrlen") = yrlen_vec,
                      Rcpp::Named("tx_cost") = tx_cost_vec,
                      Rcpp::Named("hosp_cost") = hosp_cost_vec, 
                      Rcpp::Named("mgmt_cost") = mgmt_cost_vec,
                      Rcpp::Named("si_cost") = si_cost_vec,
                      Rcpp::Named("prod_loss") = prod_loss_vec);
}

// Sample from Hernandez Alva (2013) Mixture Model
// Note the class 4 is the reference category in the paper 
//' @export
// [[Rcpp::export]]
double sim_utility_mixture1C(arma::rowvec beta1, arma::rowvec beta2, 
                arma::rowvec beta3, arma::rowvec beta4,
                double alpha1, double alpha2,
                double alpha3, double alpha4,
                double alpha, 
                double epsilon1_sd, double &epsilon2_sd,
                double epsilon3_sd, double &epsilon4_sd,
                double mu, arma::mat delta,
                arma::rowvec w, arma::rowvec x, int male){
  arma::rowvec latclass_prob = mlogit_probC(w, delta);
  int latclass = hesim::rcat1C(latclass_prob); 
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
List sim_utility_mixtureC(std::vector<int> id, std::vector<int> sim, std::vector<double> haq,
              double pain_mean, double haq_mean,
              double pain_var, double haq_var, double painhaq_cor,
              std::vector<double> age, std::vector<double> male,
              arma::mat beta1, arma::mat beta2, 
              arma::mat beta3, arma::mat beta4, 
              std::vector<double> alpha1, std::vector<double> alpha2,
              std::vector<double> alpha3, std::vector<double> alpha4,
              std::vector<double> alpha,
              std::vector<double> epsilon1_sd, std::vector<double> epsilon2_sd,
              std::vector<double> epsilon3_sd, std::vector<double> epsilon4_sd,
              std::vector<double> mu_sd, arma::cube delta){
  int N = haq.size();
  std::vector<double> pain;
  std::vector<double> util;
  pain.reserve(N);
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
    pain.push_back(rbvcnormC(haq[i], pain_mean, haq_mean, pain_var, 
                             haq_var, painhaq_cor));
    x(0) = haq[i]; x(1) = pow(haq[i], 2); x(2) = pain[i]/100; 
    x(3) = age[i]/10000; x(4) = pow(age[i]/10000, 2);
    w(0) = 1.0; w(1) = haq[i]; w(2) = pain[i]/100; w(3) = pow(pain[i]/100, 2);
    util.push_back(sim_utility_mixture1C(beta1.row(s), beta2.row(s),
                                   beta3.row(s), beta4.row(s),
                                   alpha1[s], alpha2[s], alpha3[s], alpha4[s], alpha[s],
                                   epsilon1_sd[s], epsilon2_sd[s], epsilon3_sd[s], epsilon4_sd[s],
                                   mu, delta.slice(s), w, x, male[id[i]]));
  }
  return List::create(pain, util);
}

// Wailoo 2006 HAQ to Utility Conversion
//' @export
// [[Rcpp::export]]
std::vector<double> sim_utility_wailooC(std::vector<int> sim, std::vector<int> id,
                                        std::vector<double> age, double disease_duration,
                                        std::vector<double> haq0, std::vector<double> male, 
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

// Serious Infection Cost
//' @export
// [[Rcpp::export]]
std::vector<double> sim_si_costC(std::vector<double> si, std::vector<double> &yrlen,
                               std::vector<int> &sim, std::vector<double> cost){
  int N = si.size();
  std::vector<double> si_cost_vec;
  si_cost_vec.reserve(N);
  double si_cost = 0.0;
  for (int i = 0; i < N; ++i){
    if (si[i] == 1){
      si_cost = cost[sim[i]];
    }
    else {
      si_cost = 0;
    }
    si_cost_vec.push_back(si_cost);
  }
  return si_cost_vec;
}

// QALYs
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


