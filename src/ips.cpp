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

ItreatSwitch sim_itreat_switch(std::string itreat_switch_model, int line, int therapy, int nbt,
                  int acr, int eular, double das28, double sdai, double cdai,
                  arma::rowvec acr2das28, arma::rowvec acr2sdai,arma::rowvec acr2cdai,
                  arma::rowvec nma_das28_1, arma::rowvec nma_das28_2, arma::rowvec p){
  
  ItreatSwitch sim;
  sim.tswitch = 0;
  int da_change = 0;
  int da_cat = 0;
  
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
        da_cat = get_das28_cat(das28 + da_change);
      }
      else if (itreat_switch_model == "acr-sdai-switch"){
        da_change = acr2sdai(acr);
        da_cat = get_sdai_cat(sdai + da_change);
      }
      else if (itreat_switch_model == "acr-cdai-switch"){
        da_change = acr2cdai(acr);
        da_cat = get_cdai_cat(cdai + da_change);
      }
      sim.tswitch = R::rbinom(1, p(da_cat));
    }
  else if (itreat_switch_model == "das28-switch"){ // Treatment -> DA -> Switch
      if (line == 0){
        da_cat = get_das28_cat(das28 + nma_das28_1(therapy));
      }
      else{
        da_cat = get_das28_cat(das28 + nma_das28_2(therapy));
      }
      sim.tswitch = R::rbinom(1, p(da_cat));
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


// Treatment duration by eular response 
//' @export
// [[Rcpp::export]]
double sim_duration_eular(arma::rowvec x, arma::rowvec loc_mod, double anc1_mod,
                 arma::rowvec loc_good, double anc1_good, 
                 int eular, std::string dist, double cycle_length, double si_duration,
                 double anc2_mod = 0.0, double anc2_good = 0.0){
  double surv = 0.0;
  if (si_duration < 0){
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

// Treatment duration
//' @export
// [[Rcpp::export]]
double sim_duration(arma::rowvec x, arma::rowvec loc, double anc1, std::string dist, 
                          int tswitch, double cycle_length, double si_duration,
                          double anc2 = 0.0){
  double surv = 0.0;
  if (si_duration < 0 || tswitch == 1){
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

// Simulate HAQ score
//' @export
// [[Rcpp::export]] 
List sim_haqC(arma::mat therapies,
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
             arma::mat x_mort, arma::mat logor, 
             std::string dur_dist, arma::mat x_dur, 
             arma::mat dur_loc, arma::vec dur_anc1, arma::vec dur_anc2,
             arma::mat dur_loc_mod, arma::vec dur_anc1_mod, arma::vec dur_anc2_mod,
             arma::mat dur_loc_good, arma::vec dur_anc1_good, arma::vec dur_anc2_good,
             double cycle_length, int treat_gap, int cdmards, int nbt, 
             arma::mat si_loc, arma::mat si_anc1, arma::mat si_anc2, std::string si_dist, 
             arma::mat haqdelta_loghr, int max_months){
  
  // Constants
  int n_therapies = therapies.n_cols;
  int n_sims = logor.n_rows;
  int n_pat = x_mort.n_rows;
  int maxt = 100 * (12/cycle_length);
  
  // Declarations
  int counter = 0;
  double cage_i = 0; // continuous age
  int age_i = 0;
  arma::rowvec therapies_i = therapies.row(0);
  int therapies_ij = 0;
  
  // Information to store
  std::vector<int> sim;
  std::vector<int> id;
  std::vector<double> month_vec;
  std::vector<int> therapy_vec;
  //std::vector<std::string> therapy_name;
  std::vector<int> therapy_seq;
  std::vector<int> therapy_cycle;
  std::vector<int> death_vec;
  std::vector<double> age_vec;
  std::vector<double> ttd_vec;
  std::vector<int> acr_vec;
  std::vector<int> eular_vec;
  std::vector<double> haq_vec;
  std::vector<double> ttsi_vec;
  std::vector<int> si_vec;
  std::vector<double> yrlen_vec;
  
  // "OUTER" LOOP
  for (int s = 0; s < n_sims; ++s){
    
    // "INNER" LOOP
    // Loop over Individuals
    for (int i = 0; i < n_pat; ++i){
      double month = 6;
      cage_i = age0[i];
      age_i = int(cage_i);
      double haq = haq0[i];
      if (therapies.n_rows > 1){
        therapies_i = therapies.row(i);
      }
      int t_cdmards = 0;
      ItreatSwitch sim_s_t1;
      sim_s_t1.das28 = das28_0[i];
      sim_s_t1.sdai = sdai0[i];
      sim_s_t1.cdai = cdai0[i];
      
      // Loop over therapies
      for (int j = 0; j < n_therapies + 1; ++j){
        
        // Which therapy is being used and how long has it been used for
        if (j == n_therapies - 1){
           t_cdmards = 0; // counter for cdmards and nbt which does not reset after cdmards ends and nbt begins
        }
        
        if (j == n_therapies){
          therapies_ij = nbt;
          if (therapies_i(j-1) == cdmards){
              t_cdmards = therapy_cycle[counter - 1] + 1; // note 0 otherwise
          } 
        }
        else{
          therapies_ij = therapies_i(j);
        }
        
        // H1-H3: simulate change in HAQ during initial treatment phase
        ItreatHaq sim_h_t1 = sim_itreat_haq(itreat_haq_model, j, therapies_ij, nbt,
                                        nma_acr1.slice(therapies_ij).row(s), 
                                        nma_acr2.slice(therapies_ij).row(s), 
                                        nma_dhaq1.col(therapies_ij)(s),
                                        nma_dhaq2.col(therapies_ij)(s),
                                        acr2eular.slice(s), acr2haq.row(s), 
                                        eular2haq.row(s));
        
        // S1-S6: simulate treatment switching during first 6 months
        sim_s_t1 = sim_itreat_switch(itreat_switch_model, j, therapies_ij, nbt,
                                                  sim_h_t1.acr, sim_h_t1.eular,
                                                  sim_s_t1.das28, sim_s_t1.sdai, sim_s_t1.cdai,
                                                  acr2das28.row(s), acr2sdai.row(s), acr2cdai.row(s),
                                                  nma_das28_1.row(s), nma_das28_2.row(s), 
                                                  tswitch_da.row(s));
        
        // Treatment duration
        double si_duration_j = (-cycle_length/12 + rsurvC(as_scalar(si_loc.row(s).col(therapies_ij)), 
                                              as_scalar(si_anc1.row(s).col(therapies_ij)),
                                              si_dist,as_scalar(si_anc2.row(s).col(therapies_ij))
                                              )) * (12/cycle_length);
        double duration_j = 0;
        if (itreat_switch_model == "acr-eular-haq"){
            duration_j = sim_duration_eular(x_dur.row(i), dur_loc_mod.row(s), dur_anc1_mod(s), 
                                                 dur_loc_good.row(s), dur_anc1_good(s), 
                                                 sim_h_t1.eular, dur_dist, cycle_length, si_duration_j,
                                                 dur_anc2_mod(s), dur_anc2_good(s));
        }
        else {
          duration_j = sim_duration(x_dur.row(i), dur_loc.row(s), dur_anc1(s), dur_dist, 
                                    sim_s_t1.tswitch, cycle_length, si_duration_j,
                                                 dur_anc2(s));
        }
        
        // Loop over time
        for (int t = 0; t < maxt; ++t){
          
          // Stop j loop if cycle is larger than treatment duration
          if (j < n_therapies && t > 0 && t >= duration_j + 1 + treat_gap){
            cage_i = cage_i - cycle_length/12 + 0.5;
            month = month - cycle_length + 6;
            age_i = int(cage_i);
            break;
          }
          
          // Draw death indicator
          double death = sample_deathC(age_i, male[i], lifetable_male, lifetable_female, 
                                       x_mort.row(i), logor.row(s), haq0[i], haq,
                                      cycle_length, month, haqdelta_loghr.row(s));
          
          
          // Update HAQ score
          if (t == 0){ // initial haq change
            haq = haq + sim_h_t1.dhaq; 
          }
          else if (t > 0 && t <= duration_j){
            if(cdmards_haq_model == "lcgm" && (therapies_ij == nbt | therapies_ij == cdmards)){
                haq = haq + sim_dhaq_lcgm1C(t_cdmards * cycle_length/12, cycle_length, cage_i, 1 - male[i],
                                           das28_0[i], haq_lcgm_delta.slice(s), haq_lcgm_beta.slice(s));
            } else{
              update_haq_t(haq, as_scalar(haq_lprog_therapy.row(s).col(therapies_ij)),
                           haq_lprog_age.row(s), cage_i, cycle_length);
            }
          }
          else if (t > duration_j && t < duration_j + 1 && j < n_therapies){ // rebound
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
          
          // Fill vectors
          sim.push_back(s);
          id.push_back(i);
          month_vec.push_back(month);
          therapy_vec.push_back(therapies_ij);
          //therapy_name.push_back(therapy_names[j]);
          therapy_seq.push_back(j); // should be j
          therapy_cycle.push_back(t);
          death_vec.push_back(death);
          age_vec.push_back(cage_i);
          if (therapies_ij != nbt){
            ttd_vec.push_back(duration_j - t);
          }
          else {
            ttd_vec.push_back(NA_REAL);
          }
          
          acr_vec.push_back(sim_h_t1.acr);
          eular_vec.push_back(sim_h_t1.eular);
          haq_vec.push_back(haq);
          if (therapies_ij != nbt){
            ttsi_vec.push_back(si_duration_j - t);
          }
          else{
            ttsi_vec.push_back(NA_REAL);
          }
          if (si_duration_j < duration_j && t > duration_j && t < duration_j + 1
                && therapies_ij != nbt) {
            si_vec.push_back(1);
          } 
          else if (si_duration_j < 0 && therapies_ij != nbt) {
            si_vec.push_back(1);
          }
          else {
            si_vec.push_back(0);
          }
          if (t == 0) {
            yrlen_vec.push_back(0.5);
          } 
          else {
            yrlen_vec.push_back(cycle_length/12);
          }
          
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
  return List::create(sim, id, month_vec, therapy_vec, therapy_seq, therapy_cycle, death_vec, 
                      age_vec, ttd_vec, acr_vec, eular_vec,
                      haq_vec, ttsi_vec, si_vec, yrlen_vec);
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
                                        std::vector<double> age, std::vector<double> disease_duration,
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
      b_disease_duration[sim[i]] * disease_duration[id[i]] + b_haq0[sim[i]] * haq0[id[i]] +
      b_male[sim[i]] * male[id[i]] + b_prev_dmards[sim[i]] * prev_dmards[id[i]] + 
      b_haq[sim[i]] * haq[i];
    logit_xb.push_back(1/(1 + exp(-xb)));
  }
  return logit_xb;
}

// Treatment Costs
//' @export
// [[Rcpp::export]]
List treat_costC(std::vector<int> therapy, std::vector<int> therapy_cycle,
              std::vector<double> id, std::vector<double> weight, double cycle_length,
             std::vector<double> ann_infusion_cost, std::vector<double> ann_rx_cost,
             std::vector<double> init_infusion_cost, std::vector<double> init_rx_cost,
             std::vector<int> weight_based,
             std::vector<double> ann_wgt_slope, std::vector<double> init_wgt_slope,
             std::vector<double> ann_util, std::vector<double> init_util,
             std::vector<double> strength, std::vector<double> price,
             int cdmards, int tcz, int tczmtx){
  int N = therapy.size();
  std::vector<double> treat_infusion_cost;
  std::vector<double> treat_rx_cost;
  treat_infusion_cost.reserve(N);
  treat_rx_cost.reserve(N);
  double init_wcost = 0.0;
  double ann_wcost = 0.0;
  double init_infusion_cost_j = 0.0;
  double init_rx_cost_j = 0.0;
  double ann_infusion_cost_j = 0.0;
  double ann_rx_cost_j = 0.0;
  
  // Replace annual cost with costs per cycle
  for (int k = 0; k < ann_infusion_cost.size(); ++k){
    ann_infusion_cost[k] = ann_infusion_cost[k] * cycle_length/12;
    ann_rx_cost[k] = ann_rx_cost[k] * cycle_length/12;
    ann_util[k] = ann_util[k] * cycle_length/12;
  }
  
  // Loop over simulated data
  for (int i = 0; i < N; ++i){
    int j = therapy[i];
    
    // initial 6-month cost
    if (therapy_cycle[i] == 0){
      // weight based cost
      if (weight_based[j] == 1){
        int init_wunit1 = ceil(weight[id[i]]* init_wgt_slope[j] / strength[j]);
        init_wcost = price[j] * init_util[j] * init_wunit1;
      }
      else{
        init_wcost = 0;
      }
      
      // non-weight based cost
      if (j == tcz){
        if (weight[id[i]] >= 100){
          init_infusion_cost_j = init_infusion_cost[tcz] * 2;
          init_rx_cost_j = init_rx_cost[tcz] * 2;
        }
        else {
          init_infusion_cost_j = init_infusion_cost[j];
          init_rx_cost_j = init_rx_cost[j];
        }
      }
      else if (j == tczmtx){
        if (weight[id[i]] >= 100){
          init_infusion_cost_j = init_infusion_cost[tcz] * 2 + init_infusion_cost[cdmards];
          init_rx_cost_j = init_rx_cost[tcz] * 2 + init_infusion_cost[cdmards];
        }
        else {
          init_infusion_cost_j = init_infusion_cost[j];
          init_rx_cost_j = init_rx_cost[j];
        }
      }
      else {
        init_infusion_cost_j = init_infusion_cost[j];
        init_rx_cost_j = init_rx_cost[j];
      }
      double init_rx_cost = init_rx_cost_j + init_wcost;
      treat_infusion_cost.push_back(init_infusion_cost_j);
      treat_rx_cost.push_back(init_rx_cost);
    }
    
    // Annual post 6-month cost
    else {
      // weight based cost
      if (weight_based[j] == 1){
        int ann_wunit1 = ceil(weight[id[i]]* ann_wgt_slope[j] / strength[j]);
        ann_wcost = price[j] * ann_util[j] * ann_wunit1;
      }
      else{
        ann_wcost = 0;
      }
      
      // non-weight based cost
      if (j == tcz){
        if (weight[id[i]] >= 100){
          ann_infusion_cost_j = ann_infusion_cost[tcz] * 2;
          ann_rx_cost_j = ann_rx_cost[tcz] * 2;
        }
        else {
          ann_infusion_cost_j = ann_infusion_cost[j];
          ann_rx_cost_j = ann_rx_cost[j];
        }
      }
      else if (j == tczmtx){
        if (weight[id[i]] >= 100){
          ann_infusion_cost_j = ann_infusion_cost[tcz] * 2 + ann_infusion_cost[cdmards];
          ann_rx_cost_j = ann_rx_cost[tcz] * 2 + ann_rx_cost[cdmards];
        }
        else {
          ann_infusion_cost_j = ann_infusion_cost[j];
          ann_rx_cost_j = ann_rx_cost[j];
        }
      }
      else {
        ann_infusion_cost_j = ann_infusion_cost[j];
        ann_rx_cost_j = ann_rx_cost[j];
      }      
      double ann_rx_cost =  ann_rx_cost_j  + ann_wcost;
      treat_infusion_cost.push_back(ann_infusion_cost_j);
      treat_rx_cost.push_back(ann_rx_cost);
    }
  }
  return List::create(treat_infusion_cost, treat_rx_cost);
}

// Hospitalization costs
//' @export
// [[Rcpp::export]]
std::vector<double> hosp_costC(std::vector<double> &haq, std::vector<double> &yrlen,
                               std::vector<int> &sim, arma::mat &unit_cost_mat, 
                               arma::mat &hosp_days_mat){
  int N = haq.size();
  std::vector<double> hosp_cost_vec;
  hosp_cost_vec.reserve(N);
  double hosp_days = 0.0;
  double hosp_cost = 0.0;
  for (int i = 0; i < N; ++i){
    if (haq[i] < 0.5){
      hosp_days = hosp_days_mat.row(sim[i])(0);
      hosp_cost = hosp_days * unit_cost_mat.row(sim[i])(0) * yrlen[i];
    }
    else if (haq[i] >= 0.5 && haq[i] < 1){
      hosp_days = hosp_days_mat.row(sim[i])(1);
      hosp_cost = hosp_days * unit_cost_mat.row(sim[i])(1) * yrlen[i];
    }
    else if (haq[i] >= 1 && haq[i] < 1.5){
      hosp_days = hosp_days_mat.row(sim[i])(2);
      hosp_cost = hosp_days * unit_cost_mat.row(sim[i])(2) * yrlen[i];
    }
    else if (haq[i] >= 1.5 && haq[i] < 2){
      hosp_days = hosp_days_mat.row(sim[i])(3);
      hosp_cost = hosp_days * unit_cost_mat.row(sim[i])(3) * yrlen[i];
    }
    else if (haq[i] >= 2 && haq[i] < 2.5){
      hosp_days = hosp_days_mat.row(sim[i])(4);
      hosp_cost = hosp_days * unit_cost_mat.row(sim[i])(4) * yrlen[i];
    }
    else if (haq[i] >= 2.5){
      hosp_days = hosp_days_mat.row(sim[i])(5);
      hosp_cost = hosp_days * unit_cost_mat.row(sim[i])(5) * yrlen[i];
    }
    hosp_cost_vec.push_back(hosp_cost);
  }
  return hosp_cost_vec;
}

// Productivity Loss
//' @export
// [[Rcpp::export]]
std::vector<double> prod_lossC(std::vector<double> &haq, std::vector<double> &yrlen,
                               std::vector<int> &sim, std::vector<double> beta){
  int N = haq.size();
  std::vector<double> prod_loss_vec;
  prod_loss_vec.reserve(N);
  for (int i = 0; i < N; ++i){
    prod_loss_vec.push_back(haq[i] * beta[sim[i]] * yrlen[i]);
  }
  return prod_loss_vec;
}

// Serious Infection Cost
//' @export
// [[Rcpp::export]]
std::vector<double> si_costC(std::vector<double> si, std::vector<double> &yrlen,
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

// General management cost
//' @export
// [[Rcpp::export]]
std::vector<double> mgmt_costC(std::vector<double> &yrlen,
                               std::vector<int> &sim, std::vector<double> &cost){
  int N = sim.size();
  std::vector<double> mgmt_cost_vec;
  mgmt_cost_vec.reserve(N);
  for (int i = 0; i < N; ++i){
    mgmt_cost_vec.push_back(cost[sim[i]] * yrlen[i]);
  }
  return mgmt_cost_vec;
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

