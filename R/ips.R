#' The IVI-RA simulation
#'
#' Run the IVI-RA individual patient simulation model.
#' 
#' @param arms Name of arms in the treatment treatment sequence. May be a vector consisting of 
#' a single treatment sequence or a matrix of unique sequences for each patient.
#' @param input_data An object of class 'input_data' returned from \link{get_input_data}.
#' @param pars List of sampled parameter values generated from \link{sample_pars}.
#' @param model_structure Object of class model structure generated from \link{sample_pars}.
#' @param max_months Maximum number of months to run the model for. Default is NULL which implies that
#' the model is simulated over each patient's lifetime.
#' @param treatments_lookup Vector of names of all treatments included in the parameter
#' estimates. Index of of treatments in \code{arms} are matched against treatments in
#' \code{treatments_lookup} by name. Indices of treatment-specific parameter estimates must be 
#' in the same order as treatments in \code{treatments_lookup}.   
#' @param cdmards_ind Index for cDMARDs.
#' @param nbt_ind Index for the non-biologic (i.e. a term used to define a selection of treatments clinicians use
#' after the last biologic in a treatment sequence).
#' @export 
sim_iviRA <- function(arms, input_data, pars, model_structure, 
                      max_months = NULL, treatment_lookup = iviRA::treatments$sname,
                    cdmards_ind = which(iviRA::treatments$sname == "cdmards"),
                    nbt_ind = which(iviRA::treatments$sname == "nbt"),
                    check = TRUE){
  
  # PREPPING FOR THE SIMULATION
  ## check correct object types used as arguments
  if (!inherits(input_data, "input_data")){
    stop("The argument 'input_data' must be of class 'input_data'")
  }
  if (!inherits(model_structure, "model_structure")){
      stop("The argument 'model_structure' must be of class 'model_structure'")
  }
  
  ## treatment arm indices
  if (is.vector(arms)) arms <- matrix(arms, nrow = 1)
  arminds <- matrix(match(arms, treatment_lookup), nrow = nrow(arms),
                    ncol = ncol(arms))
  
  ## default internal values
  treat_gap <- 0
  cycle_length <- 6
  if (is.null(max_months)){
    max_months <- 12 * 150
  }
  prob.switch.da <- matrix(rep(c(0, 0, 1, 1), each = pars$n), ncol = 4)
  
  ## indexing
  cdmards.ind <- cdmards_ind - 1
  nbt.ind <- nbt_ind - 1
  arminds <- arminds - 1

  ## survival parameters
  ttd.dist <- model_structure[["ttd_dist"]]
  pars.ttd.all <- pars$ttd.all[[ttd.dist]]
  pars.ttd.da <- pars$ttd.da[[ttd.dist]]
  pars.ttd.em <- pars$ttd.eular$moderate[[ttd.dist]]
  pars.ttd.eg <- pars$ttd.eular$good[[ttd.dist]]
  si.dist <- "exp"
  pars.si <- pars$ttsi
  
  # treatment costs
  tc <- pars$treat.cost
  tc.arms <- cbind(arms, "nbt")
  lookup.inds <- match(tc.arms, tc$lookup$sname)
  agents <- aperm(array(match(unlist(tc$lookup[lookup.inds, -1, with = FALSE]),
                         iviRA::treat.cost$cost$sname) - 1,
                   dim = c(nrow(tc.arms), ncol(tc.arms), ncol(tc$lookup) - 1)),
                  perm = c(2, 3, 1))
  
  # RUNNING THE SIMULATION
  sim.out <- sim_iviRA_C(arm_inds = arminds, 
                         haq0 = input_data$haq0, das28_0 = input_data$das28,
                     sdai0 = input_data$sdai, cdai0 = input_data$cdai,
                     age0 = input_data$age, male = input_data$male, 
                     prev_dmards = input_data$prev.dmards,
                     itreat_haq_model = model_structure["itreat_haq"], 
                     itreat_switch_model = model_structure["itreat_switch"],
                     nma_acr1 = pars$acr$p1, nma_acr2 = pars$acr$p2, 
                     nma_dhaq1 = pars$haq$dy1, nma_dhaq2 = pars$haq$dy2,
                     nma_das28_1 = pars$das28$dy1, nma_das28_2 = pars$das28$dy2,
                     acr2eular = pars$acr2eular, acr2haq = pars$acr2haq, eular2haq = pars$eular2haq,
                     acr2das28 = pars$acr2das28, acr2sdai = pars$acr2sdai, acr2cdai = pars$acr2cdai,
                     tswitch_da = prob.switch.da,
                     haq_lprog_therapy = pars$haq.lprog.tx, haq_lprog_age = pars$haq.lprog.age,
                     haq_lcgm_delta = pars$haq.lcgm$delta, haq_lcgm_beta = pars$haq.lcgm$beta, 
                     cdmards_haq_model = model_structure["cdmards_haq_model"],
                     rebound_factor = pars$rebound, 
                     lifetable_male = pars$lt$male, lifetable_female = pars$lt$female,
                     x_mort = input_data$x.mort, logor_mort = pars$logor.mort, 
                     ttd_dist = ttd.dist, x_ttd = input_data$x.ttd,
                     ttd_all_loc = pars.ttd.all$sample[, pars.ttd.all$loc.index, drop = FALSE],
                     ttd_all_anc1 = pars.ttd.all$sample[, pars.ttd.all$anc1.index, drop = FALSE],
                     ttd_all_anc2 = pars.ttd.all$sample[, pars.ttd.all$anc2.index, drop = FALSE],
                     ttd_da_loc = pars.ttd.da$sample[, pars.ttd.da$loc.index, drop = FALSE],
                     ttd_da_anc1 = pars.ttd.da$sample[, pars.ttd.da$anc1.index, drop = FALSE],
                     ttd_da_anc2 = pars.ttd.da$sample[, pars.ttd.da$anc2.index, drop = FALSE],
                     ttd_eular_loc_mod = pars.ttd.em$sample[, pars.ttd.em$loc.index, drop = FALSE], 
                     ttd_eular_anc1_mod = pars.ttd.em$sample[, pars.ttd.em$anc1.index, drop = FALSE],
                     ttd_eular_anc2_mod = pars.ttd.em$sample[, pars.ttd.em$anc2.index, drop = FALSE], 
                     ttd_eular_loc_good = pars.ttd.eg$sample[, pars.ttd.eg$loc.index, drop = FALSE], 
                     ttd_eular_anc1_good = pars.ttd.eg$sample[, pars.ttd.eg$anc1.index, drop = FALSE],
                     ttd_eular_anc2_good = pars.ttd.eg$sample[, pars.ttd.eg$anc2.index, drop = FALSE],
                     cycle_length = cycle_length, treat_gap = treat_gap, 
                     cdmards = cdmards.ind, nbt = nbt.ind,
                     si_loc = pars.si, 
                     si_anc1 = matrix(NA, nrow = pars$n, ncol = ncol(pars.si)), 
                     si_anc2 = matrix(NA, nrow = pars$n, ncol = ncol(pars.si)), 
                     si_dist = si.dist, 
                     haqdelta_loghr = pars$mort.loghr.haqdif, max_months = max_months,
                     hosp_days = pars$hosp.cost$hosp.days, cost_pday = pars$hosp.cost$cost.pday,
                     mgmt_cost = rowSums(pars$mgmt.cost), 
                     si_cost = pars$si.cost, prod_loss = pars$prod.loss,
                     tc_agents_ind = agents, tc = tc[c("cost", "discount")], 
                     weight = input_data$weight)
  sim.out <- as.data.table(sim.out)
  
  ## C++ to R indices
  sim.out[, sim := sim + 1]
  sim.out[, id := id + 1]
  sim.out[, tx_seq := tx_seq + 1]
  sim.out[, tx_cycle := tx_cycle + 1]
  sim.out[, tx := tx + 1]

  ## simulate utility and calculate QALYs
  if (model_structure["utility_model"] == "mixture"){
      sim.out <- cbind(sim.out, sim_utility_mixture(sim.out, male = input_data$male,
                                         pars = pars$utility.mixture))
      
  } else if (model_structure["utility_model"] == "wailoo"){
      sim.out[, utility := sim_utility_wailoo(sim.out, haq0 = input_data$haq0,
                                        male = input_data$male,
                                        prev_dmards = input_data$prev.dmards,
                                        coefs = pars$utility.wailoo)]
  }
  sim.out[, qalys := sim_qalys(simhaq = sim.out, utility = sim.out$utility, 
                               si_ul = pars$si.ul)]
  
  # RETURN
  return(sim.out)
}

#' Select model structure
#'
#' Select the model structure to be used in the IVI-RA individual patient simulation.
#' 
#' @param itreat_haq Model structure relating treatment to HAQ during the first 6 months of 
#' treatment.
#' @param itreat_switch Model structure relating treatment to switching during the first
#' 6 months of treatment.
#' @param cdmards_haq_model Model used for long-term HAQ progression.
#' @param ttd_dist Distribution used to model time to treatment discontinuaton. Options are the
#'  exponential (\code{exponential}), Weibull (\code{weibull}), Gompertz (\code{gompertz}),
#' gamma (\code{gamma}), log-logistic (\code{llogis}), lognormal (\code{lnorm}), 
#' and generalized gamma (\code{gengamma}) distributions.
#' @param utility_model Model used to estimate patient utility as a function of HAQ and patient
#' characteristics.
#' @export
select_model_structure <- function(itreat_haq = c("acr-haq", "acr-eular-haq", "haq"),
                                   itreat_switch = c("acr-switch", "acr-das28-switch",
                                                     "acr-sdai-switch", "acr-cdai-switch", 
                                                     "das28-switch", "acr-eular-switch"),
                                   cdmards_haq_model = c("lcgm", "linear"), 
                                   ttd_dist = c("exponential", "weibull", "gompertz", 
                                           "gamma", "llogis", "lnorm", "gengamma"),
                                   utility_model = c("mixture", "wailoo")){
  itreat.haq <- match.arg(itreat_haq)
  itreat.switch <- match.arg(itreat_switch)
  cdmards.haq_model <- match.arg(cdmards_haq_model)
  ttd.dist <- match.arg(ttd_dist) 
  utility.model <- match.arg(utility_model)
  if (itreat.haq == "acr-haq"){
      if (itreat.switch == "acr-eular-switch"){
        stop ("'acr-eular-switch' cannot be used when itreat_haq equals 'acr-haq'.")
      }
  }
  if (itreat.haq == "haq"){
      if (itreat.switch != "das28-switch"){
        stop("Only 'das28-swtich' can be used when itreat_haq equals 'haq'")
      }
  }
  model.structure <- c("itreat_haq" = itreat.haq, "itreat_switch" = itreat.switch,
                       "cdmards_haq_model" = cdmards.haq_model,
                       "ttd_dist" = ttd.dist,
                       "utility_model" = utility.model)
  class(model.structure) <- "model_structure"
  return(model.structure)
}

#' Simulate HAQ over time
#'
#' An individual patient simulation of HAQ scores for rheumatoid arthritis patients given a sequence of J 
#' treatments. An "outer" loop iterates over S draws from sampled parameter sets and an "inner" loop iterates of
#' N patients. Cycle lengths are 6 months long. The simulation is written in C++ for speed. 
#' 
#' @param arminds Indices of treatment arms consisting of sequences of therapies (each element
#' is the index of a treatment in \code{iviRA::treatments}. May be a vector consisting of a single treatment sequence or a matrix 
#' of unique sequences for each patient.
#' @param input_data List of input data. Required inputs are \code{haq0}, \code{age}, \code{male}, \code{x.mort}, 
#' and \code{x.ttd} as generated from \link{input_data}.
#' @param pars List of parameters. Required parameters are \code{rebound}, \code{acr1}, \code{acr2}, \code{acr2eular}, \code{eular2haq}, 
#' \code{haq.lprog.tx}, \code{haq.lprog.age}, \code{logor.mort}, \code{mort.loghr.haqdif}, \code{ttsi},
#' \code{ttd.eular.mod}, \code{ttd.eular.good}, and \code{lt} as generated from \link{sample_pars}. Additionally, if \code{cdmards_prog} is equal 
#' to "lcgm", then \code{haq.lcgm} must be included. 
#' @param itreat_haq How should the relationship between treatment and HAQ during the initial 
#' treatment phase (i.e., first 6 months) be modeled? Options, which are equivalent to H1-H3
#' in the documentation are
#' H1: Treatment -> ACR -> HAQ (\code{acr-haq}), 
#' H2: Treatment -> ACR -> EULAR -> HAQ (\code{acr-eular-haq}), or
#' H3:Treatment -> HAQ (\code{haq}). 
#' @param itreat_switch How should the relationship between treatment and switching to a new 
#' treatment during the initial treatment phase (i.e., first 6 months) be modeled. Options, which
#' are equivalent to S1-S6 in the documentation are
#' S1: Treatment -> ACR -> Switch (\code{acr-switch}), 
#' S2: Treatment -> ACR -> DAS28 -> Switch (\code{acr-das28-switch}),
#' S3: Treatment -> ACR -> SDAI -> Switch (\code{acr-sdai-switch}),
#' S4: Treatment -> ACR -> CDAI -> Switch (\code{acr-cdai-switch}),
#' S5: Treatment -> DAS28 -> Switch (\code{das28-switch}),
#' S6: Treatment -> ACR -> EULAR -> Switch (\code{acr-eular-switch}). 
#' @param haq_prog_model Model used to simulate the progression of HAQ. Options are "lcgm" and "linear", with
#' "lcgm" as the default. If "lcgm" is chosen, then a latent class growth model is used for cDMARDs
#' and NBT but a constant annual rate is is assumed for all other therapies; otherwise 
#' a constant linear HAQ progression is assumed for all therapies including cDMARDs and NBT.
#' @param ttd_dist Survival distribution for treatment duration.
#' @param si_dist Survival distribution for serious infections. Currently rate is assumed
#' to be constant so an exponential distribution is used.
#' @param max_months Maximum number of months to run the model for. Default is NULL which implies that
#' the model is simulated over each patient's lifetime.
#' @param cdmards_ind Index for cDMARDs.
#' @param nbt_ind Index for the non-biologic (i.e. a term used to define a selection of treatments clinicians use
#' after the last biologic in a treatment sequence).
#' @param check Should the function check parameters and input data passed to the model? Default is TRUE.
#' 
#' @return A \code{data.table} with the following columns:
#'\describe{
#' \item{sim}{Simulation number. Indexed from 1 to S where S is the number of randomly sampled parameter sets (e.g. n from \link{sample_pars}).}
#' \item{id}{ID number. Indexed from 1 to N where N is the number of simulated patients (e.g. n from \link{sample_pats}).}
#' \item{month}{Month since a simulated patient began the first treatment in a treatment sequence.}
#' \item{tx}{Treatment used. Given J total therapies, the first J - 1 therapies match the indices from \code{arminds}. The final \code{tx}
#'  is always the non-biologic treatment (\code{nbt_ind}).}
#' \item{tx_seq}{Number of treatment in a treatment sequence. First treatment equal to 1, 
#' second treatment equal to 2, \ldots}
#' \item{tx_cycle}{Number of model cycles since a patient began taking a given treatment in a treatment sequence. \code{tx_cycle} = 1
#' during the initial 6-month treatment period.}
#' \item{death}{Equal to 1 if patient died during the model cycle and 0 otherwise. }
#' \item{age}{Age of patient which increases with the model cycles.}
#' \item{ttd}{Time to treatment discontinuation. Measured in terms of model cycles (e.g. ttd = 2 if treatment will discontinue in 
#' 1 year given 6-months cycles). \code{ttd} is measured at the end of each cycle. Patients switch treatments during the cycle in which \code{ttd} 
#' becomes negative and HAQ rebounds during the cycle.}
#' \item{acr}{Simulated ACR response during the initial 6-month period for a new treatment Constant within \code{tx}.
#' Categories are 0 (ACR < 20), 1 (ACR 20-50), 2 (ACR 50-70), and 3 (ACR 70+).}
#' \item{eular}{Simulated EULAR response during the initial 6-month period for a new treatment Constant within \code{tx}. 
#' Categories are 0 (no EULAR response), 1 (moderate EULAR response), and 2 (good EULAR response).}
#' \item{haq}{HAQ score. Restricted to range between 0 and 3.}
#' \item{ttsi}{Time to serious infection. Like \code{ttd}, measured in terms of model cycles. \code{ttsi} is measured at the end of each 
#' cycle.}
#' \item{si}{Equal to 1 if treatment discontinuation was caused by a serious infection and 0 otherwise.}
#' \item{yrlen}{Length of a model cycle in years. Equal to 0.5 given 6-month cycles.}
#'}
#' @export
sim_haq <- function(arminds, haq0, das28_0, sdai0, cdai0, age0, male, prev_dmards,
                    itreat_haq, itreat_switch, 
                    nma_acr1, nma_acr2, nma_dhaq1, nma_dhaq2, nma_das28_1, nma_das28_2,
                    acr2eular, acr2haq, eular2haq, acr2das28, acr2sdai, acr2cdai,
                    tswitch_da, 
                    haq_lprog_therapy, haq_lprog_age,
                    haq_lcgm_delta, haq_lcgm_beta, cmdards_haq_model,
                    rebound_factor,
                    lifetable_male, lifetable_female, 
                    x_mort, logor,
                    dur_dist, x_dur,
                    ttd_all_loc, ttd_all_anc1, ttd_all_anc2,
                    ttd_da_loc, ttd_da_nac1, ttd_da_anc2,
                    ttd_eular_loc_mod, ttd_eular_anc1_mod, ttd_eular_anc2_mod,
                    ttd_eular_loc_good, ttd_eular_anc1_good, ttd_eular_anc2_good,
                    cycle_length, treat_gap, cdmards, nbt, 
                    si_loc, si_anc1, si_anc2, si_dist, haqdelta_loghr, max_months){
  
  # run simulation
  browser()
  simout <- sim_iviRA_C(arm_inds = arminds, haq0 = haq0, das28_0 = das28_0, 
                        sdai0 = sdai0, cdai0 = cdai0, age0 = age0, male = male,
                        prev_dmards = prev_dmards,
                      itreat_haq_model = itreat_haq, itreat_switch_model = itreat_switch, 
                      nma_acr1 = nma_acr1, nma_acr2 = nma_acr2, 
                      nma_dhaq1 = nma_dhaq1, nma_dhaq2 = nma_dhaq2, 
                      nma_das28_1 = nma_das28_1, nma_das28_2 = nma_das28_2,
                     acr2eular = acr2eular, acr2haq = acr2haq, eular2haq = eular2haq, 
                     acr2das28 = acr2das28, acr2sdai = acr2sdai, acr2cdai = acr2cdai,
                     tswitch_da = tswitch_da, 
                     haq_lprog_therapy = haq_lprog_therapy, haq_lprog_age = haq_lprog_age,
                     haq_lcgm_delta = haq_lcgm_delta, haq_lcgm_beta = haq_lcgm_beta, 
                     cdmards_haq_model = cmdards_haq_model,
                     rebound_factor = rebound_factor,
                     lifetable_male = lifetable_male, lifetable_female = lifetable_female, 
                     x_mort = x_mort, logor = logor,
                     dur_dist = dur_dist, x_dur = x_dur,
                     ttd_all_loc, ttd_all_anc1, ttd_all_anc2,
                     ttd_da_loc, ttd_da_nac1, ttd_da_anc2,
                     ttd_eular_loc_mod, ttd_eular_anc1_mod, ttd_eular_anc2_mod,
                     ttd_eular_loc_good, ttd_eular_anc1_good, ttd_eular_anc2_good,
                     cycle_length, treat_gap, cdmards, nbt, 
                     si_loc, si_anc1, si_anc2, si_dist, haqdelta_loghr, max_months)
  simout <- as.data.table(simout)
  colnames(simout) <- c("sim", "id", "month", "tx", "tx_seq", 
                     "tx_cycle", "death", "age", "ttd", "acr", "eular", 
                     "haq", "ttsi", "si", "yrlen")
  
  # C++ to R indices
  simout[, sim := sim + 1]
  simout[, id := id + 1]
  simout[, tx_seq := tx_seq + 1]
  simout[, tx_cycle := tx_cycle + 1]
  simout[, tx := tx + 1]
  return(simout)
}

#' Check parameters of sim_haq
#'
#' Error messages when incorrect inputs are passed to sim_haq.
#' 
#' @param input_data \code{input_data} as passed to \link{sim_haq}.
#' @param pars \code{pars} as passed to \link{sim_haq}.
#' @param itreat_haq Treatment to HAQ pathway first 6 months.
#' @param itreat_switch Treatment switching pathway first 6 months.
#' @keywords internal
check_sim_haq <- function(arminds, input_data, pars, itreat_haq, itreat_switch){
  names.dist <- c("exponential", "exp", "weibull", "gompertz", "gamma", "llogis",
                  "lnorm", "gengamma")
  arminds.unique <- unique(c(arminds))
  
  # check input data
  if(is.null(input_data$haq0)) stop("'haq0' element of input_data list not given")
  if(is.null(input_data$age)) stop("'age' element of input_data list not given")
  if(is.null(input_data$male)) stop("'male' element of input_data list not given")
  if(is.null(input_data$x.mort)) stop("'x.mort' element of input_data list not given")
  if(is.null(input_data$x.ttd)) stop("'x.ttd' element of input_data list not given")
  n <- input_data$n
  if(length(input_data$age) != n  |
     nrow(input_data$x.mort) !=n | nrow(input_data$x.ttd) != n) {
      stop(paste0("Number of patients not consistent accross elements of the input_data list.",
                  " Should equal ", n))
  }
  
  # check parameter inputs
  n <- pars$n
  
  ## rebound
  if(is.null(pars$rebound)) stop("'rebound' element of pars list not given")
  if(!is.vector(pars$rebound))  stop("'rebound' element of pars list must be a vector")
  if(length(pars$rebound) != n) {
    stop(paste0("Length of 'rebound' element of pars must be equal to the number of ",
         "sampled parameter sets which equals ", n))
  }
  
  ## acr2eular
  if(is.null(pars$acr2eular)) stop("'acr2eular' element of pars list not given")
  if(!is.array(pars$acr2eular))  stop("'acr2eular' element of pars list must be an array")
  if(dim(pars$acr2eular)[1] != 4) stop("First dimension of 'acr2eular' element of pars must be 4")
  if(dim(pars$acr2eular)[2] != 3) stop("Second dimension of 'acr2eular' element of pars must be 3")
  if(dim(pars$acr2eular)[3] != n) {
    stop(paste0("Third dimension of 'acr2eular' element of pars must be equal to the number of",
         "sampled parameter sets which equals ", n))
  }

  ## eular2haq
  if(is.null(pars$eular2haq)) stop("'eular2haq' element of pars list not given")
  if(!is.matrix(pars$eular2haq))  stop("'eular2haq' element of pars list must be a matrix")
  if(ncol(pars$eular2haq) != 3) stop("Number of columns of 'eular2haq' element of pars must be 3")
  if(nrow(pars$eular2haq) != n){
    stop(paste0("Number of rows in 'eular2haq' element of pars must be equal to the number ",
         "of sampled parameter sets which equals ", n))
  } 
  
  ## haq.lprog.tx
  if(is.null(pars$haq.lprog.tx)) stop("'haq.lprog.tx' element of pars list not given")
  if(!is.matrix(pars$haq.lprog.tx))  stop("'haq.lprog.tx' element of pars list must be a matrix")
  if(nrow(pars$haq.lprog.tx) != n) {
    stop(paste0("Number of rows in 'haq.lprog.tx' element of pars must be equal to the number ",
                "of sampled parameter sets which equals ", n))
  } 
  
  ## haq.lprog.age
  if(is.null(pars$haq.lprog.age)) stop("'haq.lprog.age' element of pars list not given")
  if(!is.matrix(pars$haq.lprog.age))  stop("'haq.lprog.age' element of pars list must be a matrix")
  if(ncol(pars$haq.lprog.age) != 3) stop("Number of columns in 'haq.lprog.age' element of pars must be 3")
  if(nrow(pars$haq.lprog.age) != n) {
    stop(paste0("Number of rows in 'haq.lprog.age' element of pars must be equal to the number ",
                "of sampled parameter sets which equals ", n))
  }

  ## logor.mort
  if(is.null(pars$logor.mort)) stop("'logor.mort' element of pars list not given")
  if(!is.matrix(pars$logor.mort))  stop("'logor.mort' element of pars list must be a matrix")
  if(nrow(pars$logor.mort) != n) {
    stop(paste0("Number of rows in 'logor.mort' element of pars must be equal to the number ",
                "of sampled parameter sets which equals ", n))
  }
  
  ## mort.loghr.haqdif
  if(is.null(pars$mort.loghr.haqdif)) stop("'mort.loghr.haqdif' element of pars list not given")
  if(!is.matrix(pars$mort.loghr.haqdif))  stop("'mort.loghr.haqdif' element of pars list must be a matrix")
  if(ncol(pars$mort.loghr.haqdif) != 5) stop("Number of columns in 'mort.loghr.haqdif' element of pars must be 5")
  if(nrow(pars$mort.loghr.haqdif) != n) {
    stop(paste0("Number of rows in 'mort.loghr.haqdif' element of pars must be equal to the number ",
                "of sampled parameter sets which equals ", n))
  }
  
  ## ttsi
  if(is.null(pars$ttsi)) stop("'ttsi' element of pars list not given")
  if(!is.matrix(pars$ttsi))  stop("'ttsi' element of pars list must be a matrix")
  
  ## ttd.eular
  if(is.null(pars$ttd.eular)) stop("'ttd.eular' element of pars list not given")
  for (i in 1:length(pars$ttd.eular)){
    if(!is.list(pars$ttd.eular[[i]]))  stop("Second level of 'ttd.eular' element of pars list 
                                          must be a list")
    J <- length(pars$ttd.eular[[i]])
    # for (j in 1:J){
    #   if(all(!names(pars$ttd.eular.mod) %in% names.dist)) {
    #     dists.bad <- names(pars$ttd.eular.mod)[which(!names(pars$ttd.eular.mod) %in% names.dist)]
    #     stop(paste0("sim_haq does not support at least 1 of the survival distributions (",
    #                 dists.bad,
    #                 ") that are contained in the 'ttd.eular.mod'",
    #                 " element of pars"))
    #   }
    # }  # code similar to this needs to be added to all ttd survival distributions
  }
  ## ttd.da
  if(is.null(pars$ttd.da)) stop("'ttd.da' element of pars list not given")
  
  ## lt 
  if(is.null(pars$lt)) stop("'lt' element of pars list not given")
  if (!all(names(pars$lt) %in% c("female", "male"))) stop(paste0("'lt' element of pars must contain ",
                                                          "lifetables named 'male' and 'female'"))
  if(ncol(pars$lt$male) != 3 | ncol(pars$lt$female) != 3){
    stop("Number of columns in 'lt$female' element of pars and 'lt$male' must be 3")
  } 
  
  ## NMA ACR
  nma.acr1 <- pars$acr$p1[, , arminds.unique]
  nma.acr2 <- pars$acr$p2[, , arminds.unique]
  if (itreat_haq %in% c("acr-haq", "acr-eular-haq")){
    if (any(is.na(nma.acr1) | is.na(nma.acr2))){
      stop ("'acr' element of pars has missing parameter values for one of the selected treatment 
            arms; that is, the NMA results are missing.")
    }
  }

  ## NMA HAQ
  nma.haq1 <- pars$haq$dy1[, arminds.unique]
  nma.haq2 <- pars$haq$dy2[, arminds.unique] 
  if (itreat_haq == "haq"){
    if (any(is.na(nma.haq1) | is.na(nma.haq2))){
      stop ("'haq' element of pars has missing parameter values for one of the selected treatment 
            arms; that is, the NMA results are missing.")
    }
  }
  
  ## NMA DA28
  nma.das28.1 <- pars$das28$dy1[, arminds.unique]
  nma.das28.2 <- pars$das28$dy2[, arminds.unique] 
  if (itreat_switch == "das28-switch"){
    if (any(is.na(nma.das28.1) | is.na(nma.das28.1))){
      stop ("'das28' element of pars has missing parameter values for one of the selected treatment 
            arms; that is, the NMA results are missing.")
    }
  }
}


#' Simulate utility from simulated HAQ score (Hernandez-Alva)
#'
#' Simulate utility from HAQ score simulated using \code{sim_haq} using mixture model from
#' Hernandez-Alva (2013).
#' 
#' @param simhaq Simulation output from \link{sim_haq}.
#' @param male Vector indiciating patient gender (1 = male, 0 = female).
#' @param pars List of parameters needed to simulate utility using the Hernandez (2013) mixture 
#' model. These are the parameters 'utility.mixture' described in the documentation 
#' to \link{sample_pars}.
#' @param check Should the function check parameters and input data passed to the model? Default is TRUE.
#' 
#' @return Matrix. First column is simulated pain score and second column is simulated utility. 
#' Each row corresponds to a unique patient and time-period (i.e. month) from \link{sim_haq}.
#'
#' @export
sim_utility_mixture <- function(simhaq, male, pars, check = TRUE){
  if (check) check_sim_utility_mixture(simhaq, male, pars)
  util <- sim_utility_mixtureC(simhaq$id - 1, simhaq$sim - 1, simhaq$haq, 
                               pars$pain$pain.mean, pars$pain$haq.mean, 
                               pars$pain$pain.var, pars$pain$haq.var, 
                               pars$pain$painhaq.cor, 
                               simhaq$age, male,
                               pars$beta1, pars$beta2, pars$beta3,
                               pars$beta4, 
                               pars$alpha1, pars$alpha2, pars$alpha3, 
                               pars$alpha4,
                               pars$alpha,
                               pars$epsilon1, pars$epsilon2, pars$epsilon3, 
                               pars$epsilon4,
                               pars$mu, pars$delta)
  util <- as.data.table(util)
  colnames(util) <- c("pain", "utility")
  return(util)
}

#' Check parameters of sim_utility_mixture
#'
#' Error messages when incorrect inputs are passed to sim_utility_mixture.
#' 
#' @param simhaq Simulation output from \link{sim_haq}
#' @param male Vector indiciating patient gender (1 = male, 0 = female)
#' @param pars \code{pars} as passed to \link{sim_utility_mixture}
#' @keywords internal
check_sim_utility_mixture <- function(simhaq, male, pars){
  n <- max(simhaq$sim)
  
  ## simhaq
  if(is.null(simhaq$id)) stop("'id' column of 'simhaq' is missing")
  if(is.null(simhaq$sim)) stop("'sim' column of 'simhaq' is missing")
  if(is.null(simhaq$haq)) stop("'haq' column of 'simhaq' is missing")
  
  ## male
  if(is.null(male)) stop("Object 'male' not given")
  if(!is.vector(male)) stop("'male' must be a vector")
  if(any(!male %in% c(0, 1))) stop("'male' must be equal to 0 or 1")
  
  ## pars
  # pain.mean, haq.mean, pain.var, haq.var, painhaq.cor
  for (i in c("pain.mean", "haq.mean", "pain.var", "haq.var", "painhaq.cor")){
    if(is.null(pars$pain[[i]]))  stop(paste0("'", i, "'", " element of pars list not given"))
    if(!is.vector(pars$pain$pain.mean)) stop(paste0("'", i, "'", " element of pars must be a vector"))
    if(length(pars$pain$pain.mean) > 1) stop(paste0("'", i, "'", " element of pars must be length 1"))
  }
  
  # beta1 - beta4
  for (b in c("beta1", "beta2", "beta3", "beta4")){
    if(is.null(pars[[b]])) stop(paste0("'", b, "'", " element of pars list not given"))
    if(!is.matrix(pars[[b]])) stop(paste0("'", b, "'", " element of pars list must be a matrix"))
    if(ncol(pars[[b]]) != 5) stop(paste0("Number of columns in ","'", b, "'", " element of",
                                  " pars list must be 5"))
    if(nrow(pars[[b]]) != n) {
      stop(paste0("Number of rows in ", "'", b, "'", 
                  " element of pars must be equal to the number of ",
                  "sampled parameter sets from simhaq which equals ", n, "."))
    } 
  }
  
  # alpha1 - alpha4, alpha, epsilon1 - epsilon4, mu
  for (a in c("alpha1", "alpha2", "alpha3", "alpha4", "alpha", "epsilon1",
              "epsilon2", "epsilon3", "epsilon4")){
    if(is.null(pars[[a]])) stop(paste0("'", a, "'", " element of pars list not given"))
    if(!is.vector(pars[[a]])) stop(paste0("'", a, "'", " element of pars list must be a vetor"))
    if(length(pars[[a]]) != n) {
      stop(paste0("Length of the ", "'", a, "'", 
                  " element of pars must be equal to the number of ",
                  "sampled parameter sets from simhaq which equals ", n, "."))
    } 
  }
  
  # delta
  if(is.null(pars$delta)) stop("'delta' element of pars list not given")
  if(!is.array(pars$delta)) stop("'delta' element of pars list must be an arrary")
  if(dim(pars$delta)[1] != 3) stop("First dimension of 'delta' element of pars must be 3")
  if(dim(pars$delta)[2] != 4) stop("First dimension of 'delta' element of pars must be 4")
  if(dim(pars$delta)[3] != n) {
    stop(paste0("Third dimension of 'delta' element of pars must be equal to the number ",
                "of sampled parameter sets from simhaq which equals ", n))
  }
}

#' Simulate utility from simulated HAQ score using Wailoo (2006) model
#'
#' Simulate utility from HAQ score simulated using \code{sim_haq} using model from
#' Wailoo (2006).
#' 
#' @param simhaq Simulation output from \link{sim_haq}. Variables needed are \code{sim}, \code{id},
#'  \code{age}, and \code{haq}.
#' @param haq0 HAQ score at baseline.
#' @param male Indicator = 1 for males and 0 for females.
#' @param prev_dmards Number of previous DMARDs.
#' @param coefs Matrix of coefficients needed to simulate utility using the Wailoo (2006) model. 
#' See the documentation in \link{sample_pars} for details. Note that the matrix columns must be in 
#' the same order as generated by \link{sample_pars}. 
#' @param check Should the function check parameters and input data passed to the model? Default is TRUE.
#' 
#' @details Note that disease duration is set to 18.65 years, which
#' is the mean value from the Wailoo (2006) paper used for the parameter estimates. Age and 
#' the HAQ score are taken from the simulation output.
#' @return Vector of utility scores. 
#'
#' @export
sim_utility_wailoo <- function(simhaq, haq0, male, prev_dmards,
                               coefs, check = TRUE){
  if (check) check_sim_utility_wailoo(simhaq, haq0, male, prev_dmards, coefs)
  DIS.DUR <- 18.65 # based on mean disease duration in Wailoo (2006)
  util <- sim_utility_wailooC(simhaq$sim - 1, simhaq$id - 1, simhaq$age,
                              DIS.DUR, haq0, male, 
                              prev_dmards, simhaq$haq,
                              coefs[, "int"], coefs[, "age"], coefs[, "dis_dur"], 
                              coefs[, "haq0"], coefs[, "male"],
                              coefs[, "prev_dmards"], coefs[, "haq"])
  return(util)
}

#' Check parameters of sim_utility_wailoo
#'
#' Error messages when incorrect inputs are passed to sim_utility_wailoo.
#' 
#' @param simhaq Simulation output from \link{sim_haq}
#' @param male Vector indiciating patient gender (1 = male, 0 = female)
#' @param pars \code{pars} as passed to \link{sim_utility_mixture}
#' @keywords internal
check_sim_utility_wailoo <- function(simhaq, haq0, male, prev_dmards, coefs){
  # simhaq
  if(is.null(simhaq$id)) stop("'id' column of 'simhaq' is missing")
  if(is.null(simhaq$sim)) stop("'sim' column of 'simhaq' is missing")
  if(is.null(simhaq$age)) stop("'age' column of 'simhaq' is missing")
  if(is.null(simhaq$haq)) stop("'haq' column of 'simhaq' is missing")
  
  # input_data
  if(is.null(haq0)) stop("'haq0' element of input_data list not given")
  if(is.null(male)) stop("'male' element of input_data list not given")
  if(is.null(prev_dmards)) stop("'prev.dmards' element of input_data list not given")
  n <- length(haq0)
  if(length(haq0) != n | length(male) != n | 
     length(prev_dmards) != n) {
    stop(paste0("Number of patients not consistent accross elements of the input_data list.",
                "Should equal ", n))
  }
  
  # pars
  n <- max(simhaq$sim)
  if(ncol(coefs) != 7) stop("Number of columns in 'pars' must be 7")
  if(nrow(coefs) != n) stop(paste0("Number of rows in 'pars' must be equal to the number of ",
                           "sampled parameter sets which equals ", n))
  if(any(!colnames(coefs) %in% c("int", "age", "dis_dur", "haq0", "male", "prev_dmards", "haq"))){
    stop(paste0("Column names of 'pars' must be (in order) 'int', 'age', 'dis_dur'",
                "'haq0', 'male', 'prev_dmards', and 'haq'."))
  }
}

#' Check parameters of sim_hc_cost
#'
#' Error messages when incorrect inputs are passed to sim_hc_cost.
#' 
#' @param simhaq Simulation output from \link{sim_haq}.
#' @param weight Weight of patient.
#' @param pars \code{pars} as passed to \link{sim_hc_cost}.
#' @keywords internal
check_sim_hc_cost <- function(simhaq, weight, pars){
  # simhaq
  if(is.null(simhaq$tx)) stop("'tx' column of 'simhaq' is missing")
  if(is.null(simhaq$tx_cycle)) stop("'tx_cycle' column of 'simhaq' is missing")
  if(is.null(simhaq$id)) stop("'id' column of 'simhaq' is missing")
  if(is.null(simhaq$haq)) stop("'haq' column of 'simhaq' is missing")
  if(is.null(simhaq$yrlen)) stop("'yrlen' column of 'simhaq' is missing")
  if(is.null(simhaq$sim)) stop("'sim' column of 'simhaq' is missing")
  
  ## weight
  if(is.null(weight)) stop("'weight' element not given")
  
  ## pars
  n <- pars$n
  # treatment costs
  if(is.null(pars$treat.cost)) stop("'treat.cost' element of pars not given")
  if(!is.data.frame(pars$treat.cost)) top("'treat.cost' element of pars must be a data.frame")
  for (i in c("ann_infusion_cost", "ann_rx_cost", "init_infusion_cost", "init_rx_cost",
              "weight_based", "ann_wgt_slope", "init_wgt_slope", "ann_util", 
              "init_util", "strength", "price")){
    if(is.null(pars$treat.cost[[i]])) {
      stop(paste0("'treat.cost' element of pars must contain column ", "'", i, "'"))
    } 
  }
  
  # hospital cost
  if(is.null(pars$hosp.cost)) stop("'hosp.cost' element of pars not given")
  if(!any(names(pars$hosp.cost) %in% "hosp.days")){
    stop("'hosp.days' element of pars$hosp.cost not given")
  }
  if(!is.matrix(pars$hosp.cost$hosp.days)){
    stop("'hosp.days' element of pars$hosp.cost must be a matrix")
  }
  if(ncol(pars$hosp.cost$hosp.days) != 6){
    stop("Number of columns in 'hosp.days' element of pars$hosp.cost must be 6")
  }
  if(nrow(pars$hosp.cost$hosp.days) != n){
    stop(paste0("Number of rows in 'hosp.days' element of pars$hosp.cost must ",
                "be equal to the number of sampled parameter sets which equals ", n))
  }
  if(!any(names(pars$hosp.cost) %in% "cost.pday")){
    stop("'cost.pday' element of pars$hosp.cost not given")
  }
  if(!is.matrix(pars$hosp.cost$cost.pday)){
    stop("'cost.pday' element of pars$hosp.cost must be a matrix")
  }
  if(ncol(pars$hosp.cost$cost.pday) != 6){
    stop("Number of columns in 'cost.pday' element of pars$hosp.cost must be 6")
  }
  if(nrow(pars$hosp.cost$cost.pday) != n){
    stop(paste0("Number of rows in 'cost.pday' element of pars$hosp.cost must ",
                "be equal to the number of sampled parameter sets which equals ", n))
  }
  
  # management cost
  if(is.null(pars$mgmt.cost)) stop("'mgmt.cost' element of pars not given")
  if(!is.matrix(pars$mgmt.cost)) stop("'mgmt.cost' element of pars must be a matrix")
  if(ncol(pars$mgmt.cost) != 4) stop("Number of columns in 'mgmt.cost' element of pars must be 4")
  if(nrow(pars$mgmt.cost) != n){
    stop(paste0("Number of rows of 'mgmt.cost' element of pars must be equal to the number ",
                "of sampled parameter sets which equals ", n))
  }
  
  # serious infection cost
  if(is.null(pars$si.cost)) stop("'si.cost' element of pars not given")
  if(!is.vector(pars$si.cost)) stop("'si.cost' element of pars must be a vector")
  if(length(pars$si.cost) != n) {
    stop(paste0("Length of 'si.cost' element of pars must be equal to the number ",
                "of sampled parameter sets which equals ", n))
  }
}

#' Check parameters of prod_loss
#'
#' Error messages when incorrect inputs are passed to prod_loss.
#' 
#' @param pl_haq \code{pl_haq} as passed to \link{prod_loss}
#' @keywords internal
check_prod_loss <- function(simhaq, pl_haq){
  # simhaq
  if(is.null(simhaq$haq)) stop("'haq' column of 'simhaq' is missing")
  if(is.null(simhaq$yrlen)) stop("'yrlen' column of 'simhaq' is missing")
  if(is.null(simhaq$sim)) stop("'sim' column of 'simhaq' is missing")
  
  # pl_haq
  n <- max(simhaq$sim)
  if(is.null(pl_haq)) stop("'pl_haq' is missing")
  if(!is.vector(pl_haq)) stop("'pl_haq' must be a vector")
  if(length(pl_haq) != n){
    stop(paste0("Number of rows in 'pl_haq' must ",
                "be equal to the number of sampled parameter sets which equals ", n))
  }
}

#' Simulate QALYs
#'
#' Simulate QALYs using using output produced from \code{sim_haq}.
#' 
#' @param simhaq Simulation output from \link{sim_haq}. Must include columns \code{yrlen} for
#' year length of model cycle, \code{sim} for simulation number, and \code{si} for whether a serious
#' infection occured during the model cycle. 
#' @param utility Simulated utility from \link{sim_utility_mixture} or \link{sim_utility_wailoo}.
#' @param si_ul Sampled utility loss. Equivalent to output \code{si.ul} from \link{sample_pars}.
#' @return Vector of QALYs for each simulated patient and time-period.
#'
#' @export
sim_qalys <- function(simhaq, utility, si_ul, check = TRUE){
  if (check) check_sim_qalys(simhaq, utility, si_ul)
  qalys <- qalysC(utility, simhaq$yrlen, simhaq$sim - 1, simhaq$si, si_ul)
  return(qalys)
}

#' Check parameters of sim_qalys
#'
#' Error messages when incorrect inputs are passed to sim_qalys
#' 
#' @param si_ul \code{si_ul} as passed to \link{sim_qalys}
#' @keywords internal
check_sim_qalys <- function(simhaq, utility, si_ul){
  # utility 
  if(is.null(utility)) stop("'utility' is missing")
  
  # simhaq
  if(is.null(simhaq$yrlen)) stop("'yrlen' column of 'simhaq' is missing")
  if(is.null(simhaq$sim)) stop("'sim' column of 'simhaq' is missing")
  if(is.null(simhaq$si)) stop("'si' column of 'simhaq' is missing")
  
  # si_ul
  n <- max(simhaq$sim)
  if(is.null(si_ul)) stop("'si_ul' is missing")
  if(!is.vector(si_ul)) stop("'si_ul' must be a vector")
  if(length(si_ul) != n){
    stop(paste0("Number of rows in 'si_ul' must ",
                "be equal to the number of sampled parameter sets which equals ", n))
  }
}


