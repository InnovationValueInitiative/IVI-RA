#' The IVI-RA simulation
#'
#' Run the IVI-RA individual patient simulation model.
#' 
#' @param tx_seqs Either a vector consisting of 
#' a single treatment sequence or a matrix of unique treatment sequences for each patient.
#' @param input_data An object of class 'input_data' returned from \link{get_input_data}.
#' @param pars An object of class 'par_sample' returned from \link{sample_pars}.
#' @param model_structures An object of class "model_structures" 
#' returned from \link{select_model_structures}.
#' @param max_months Maximum number of months to run the model for. Default is NULL which implies that
#' the model is simulated over each patient's lifetime.
#' @param tx_data Dataset of treatments with columns names equivalent to \code{iviRA::treatments}.
#' The indices of of treatments in \code{tx_seqs} are matched against treatments in
#' \code{tx_data$sname} by name. Indices of treatment-specific parameter estimates must be 
#' in the same order as treatments in \code{tx_data}. 
#' @param hist Is the patient tDMARD naive or tDMARD experienced?
#' @param output Specifies the format of output returned from the simulation. Options are \code{data} 
#' and \code{summary}. When \code{data} is specified, each simulated value (i.e, by model,
#' sampled parameter set, patient, and time-period) is returned in a \code{data.table}. If 
#' \code{summary} is selected, then only summary measures are returned.
#' @param discount_qalys Discount rate for QALYs. Only used when \code{output = "summary"}; 
#' otherwise, discounts can be applied to the simulated output.
#' @param discount_cost Discount rate for cost variables. Only used when \code{output = "summary"};
#' otherwise, discounts can be applied to the simulated output.
#' 
#' @return 
#' When \code{output = "data"} is selected, the simulation returns a \code{data.table} with 
#' simulated output for every model structure (\code{model}), sampled parameter set 
#' (\code{sim}), patient (\code{id}), treatment within a treatment sequence (\code{tx}),
#' and time-period (\code{month}) is returned. For more details on the output, 
#' see the 'data' section below. 
#' 
#' However, since all simulated output is returned, the size of the output can be very 
#' large and can cause memory management issues. The \code{output = "summary"}, which only
#'  provides summaries of the simulation output, can be useful in the cases. In particular, 
#'  when \code{output = "summary"}, the simulation returns three \code{data.tables}:
#' \describe{
#' \item{means}{A \code{data.table} of mean model outcomes by model structure (\code{model})
#'  and sampled parameter set (\code{sim}). For details on the variables returned see the 
#'  'means' section below.}
#' \item{time.means}{A \code{data.table} of mean model outcomes by model structure (\code{model}),
#' sampled parameter set (\code{sim}), and month (\code{month}). For details on the variables
#' returned see the 'time.means' section below.}
#' \item{out0}{Simulated model outcomes during model cycle 0 for each model structure 
#' (\code{model}), sampled parameter set (\code{sim}), patient (\code{id}), and 
#' treatment within a treatment sequence (\code{tx}). For details on the variables returned, see
#' the 'out0' section below. }
#' }
#' @section data:
#'\describe{
#' \item{model}{Integer denoting a unique model structure corresponding to the row number in
#' \code{model_structures}.}
#' \item{sim}{Simulation number denoting a randomly sampled parameter set from \link{sample_pars}.}
#' \item{id}{ID number denoting a simulated patients (e.g., from \link{sample_pop}).}
#' \item{month}{Month since a simulated patient began the first treatment in a treatment sequence.}
#' \item{tx}{Treatment used. Given J total therapies, the first J - 1 therapies match the indices from
#'  \code{tx_data}. The final \code{tx} is always the non-biologic treatment (\code{nbt_ind}).}
#' \item{line}{Line of treatment in a treatment sequence. First treatment equal to 1, 
#' second treatment equal to 2, \ldots}
#' \item{tx_cycle}{Number of model cycles since a patient began taking a given treatment in a treatment sequence. \code{tx_cycle} = 1
#' during the initial 6-month treatment period.}
#' \item{death}{Equal to 1 if patient died during the model cycle and 0 otherwise. }
#' \item{age}{Age of patient which increases with the model cycles.}
#' \item{ttd}{Time to treatment discontinuation following the initial 6 month treatment period (i.e., after
#' \code{tx_cycle=1}. Measured in terms of model cycles (e.g. ttd = 2 if treatment will discontinue in 
#'  1 year given 6-months cycles). \code{ttd} is measured at the end of each cycle.
#'  Patients switch treatments during the cycle in which \code{ttd} becomes negative and HAQ
#'   rebounds during the cycle.}
#' \item{acr}{Simulated ACR response during the initial 6-month period for a new treatment Constant within \code{tx}.
#' Categories are 0 (ACR < 20), 1 (ACR 20-50), 2 (ACR 50-70), and 3 (ACR 70+).}
#' \item{eular}{Simulated EULAR response during the initial 6-month period for a new treatment Constant within \code{tx}. 
#' Categories are 0 (no EULAR response), 1 (moderate EULAR response), and 2 (good EULAR response).}
#' \item{das28}{Disease Activity Score 28 (DAS28).}
#' \item{sdai}{Simplified Disease Activity Index (SDAI).}
#' \item{cdai}{Clinical Disease Activity Index (CDAI).}
#' \item{haq}{HAQ score. Restricted to range between 0 and 3.}
#' \item{ttsi}{Time to serious infection. Like \code{ttd}, measured in terms of model cycles. \code{ttsi} is measured at the end of each 
#' cycle.}
#' \item{si}{Equal to 1 if treatment discontinuation was caused by a serious infection and 0 otherwise.}
#' \item{yrlen}{Length of a model cycle in years. Equal to 0.5 given 6-month cycles.}
#' \item{tx_cost}{Treatment costs after discounts and rebates.}
#' \item{hosp_cost}{Hospitalization costs.}
#' \item{mgmt_cost}{General management costs.}
#' \item{si_cost}{Costs due to serious infections.}
#' \item{prod_loss}{Productivity loss (i.e., lost earnings).}
#' \item{utility}{Simulated utility score.}
#' \item{qalys}{Quality-adjusted life-years (QALYs).}
#'}
#'
#' @section means:
#' \describe{
#' \item{model}{Integer denoting a unique model structure corresponding to the row number in
#' \code{model_structures}.}
#' \item{sim}{Simulation number denoting a randomly sampled parameter set from \link{sample_pars}.}
#' \item{lys}{Life-years.}
#' \item{dlys}{Discounted life-years.}
#' \item{lys_infusion}{Total life-years with a treatment administered by infusion.}
#' \item{lys_injection}{Total life-years with a treatment administered by injection.}
#' \item{lys_oral}{Total life-years with a treatment administered orally.}
#' \item{dhaq}{Change in HAQ from baseline (i.e., month 0) to final model cycle.}
#' \item{si}{Number of serious infections.}
#' \item{qalys}{Quality-adjusted life-years (QALYs).}
#' \item{dqalys}{Discounted QALYs.}
#' \item{tx_cost}{Treatment costs.}
#' \item{dtx_cost}{Discounted treatment costs.}
#' \item{hosp_days}{Hospital days.}
#' \item{hosp_cost}{Hospital costs.}
#' \item{dhosp_cost}{Discounted hospital costs.}
#' \item{mgmt_cost}{General management costs.}
#' \item{dmgmt_cost}{Discounted general management costs.}
#' \item{si_cost}{Costs due to serious infections.}
#' \item{dsi_cost}{Discounted costs due to serious infections.}
#' \item{prod_loss}{Productivity loss (i.e., lost earnings).}
#' \item{dprod_loss}{Discounted productivity loss.}
#' \item{dhc_cost}{Discounted formal healthcare sector costs.}
#' \item{dtot_cost}{Discounted total (formal healthcare sector + productivity losses) 
#' costs.}
#' \item{yrs_since_approval}{Weighted mean of the number of years since approval of all treatments 
#' in a treatment sequence, with weights for a given treatment equal to the number of months a
#'  simulated patient used that treatment.}
#' \item{dqalys_ann}{Annualized discounted QALYs. Assumes that QALYs accrue at a constant 
#' annual rate and are calculated over the maximum number of years that each simulated patient
#' could survive during the model.}
#' \item{dhc_cost_ann}{Annualized discounted formal healthcare sector costs. Calculated in the same
#' way as \code{dqalys_ann}.}
#' \item{dprod_loss_ann}{Annualized discounted productivity losses. Calculated in the same
#' way as \code{dqalys_ann}.}
#' }
#' 
#' 
#' @section time.means:
#' \describe{
#' \item{model}{Integer denoting a unique model structure corresponding to the row number in
#' \code{model_structures}.}
#' \item{sim}{Simulation number denoting a randomly sampled parameter set from \link{sample_pars}.}
#' \item{month}{Month since a simulated patient began the first treatment in a treatment sequence.}
#' \item{alive}{Number of simulated patients alive.}
#' \item{qalys}{Quality-adjusted life-years.}
#' \item{haq}{HAQ score. Restricted to range between 0 and 3.}
#' \item{tx_cost}{Treatment costs.}
#' \item{hosp_cost}{Hospital costs.}
#' \item{mgmt_cost}{General management costs.}
#' \item{si_cost}{Costs due to serious infections.}
#' \item{prod_loss}{Productivity loss (i.e., lost earnings).}
#' }
#' 
#' @section out0:
#' \describe{
#' \item{model}{Integer denoting a unique model structure corresponding to the row number in
#' \code{model_structures}.}
#' \item{sim}{Simulation number denoting a randomly sampled parameter set from \link{sample_pars}.}
#' \item{id}{ID number denoting a simulated patients (e.g., from \link{sample_pop}).}
#' \item{tx}{Treatment used. Given J total therapies, the first J - 1 therapies match the indices from
#' }
#' \item{acr}{Simulated ACR response during the initial 6-month period for a new treatment Constant within \code{tx}.
#' Categories are 0 (ACR < 20), 1 (ACR 20-50), 2 (ACR 50-70), and 3 (ACR 70+).}
#' \item{eular}{Simulated EULAR response during the initial 6-month period for a new treatment Constant within \code{tx}. 
#' Categories are 0 (no EULAR response), 1 (moderate EULAR response), and 2 (good EULAR response).}
#' \item{ttd}{Time to treatment discontinuation. Measured in terms of model cycles (e.g. ttd = 2 if treatment will discontinue in 
#' 1 year given 6-months cycles). \code{ttd} is measured at the end of each cycle. Patients switch treatments during the cycle in which \code{ttd} 
#' becomes negative and HAQ rebounds during the cycle.}
#' \item{ttsi}{Time to serious infection. Like \code{ttd}, measured in terms of model cycles. \code{ttsi} is measured at the end of each 
#' cycle.}
#' }
#'
#' @examples 
#' pop <- sample_pop(n = 10, type = "homog")
#' tx.seq <- c("adamtx", "cdmards")
#' mod.structs <- select_model_structures(tx_ihaq = c("acr-haq", "acr-eular-haq"),
#'                                       tx_iswitch = c("acr-switch", "acr-eular-switch"),
#'                                       cdmards_haq_model = c("lcgm", "linear"),
#'                                       ttd_cause = c("all", "si"),
#'                                       ttd_dist = c("gengamma", "exponential"),
#'                                       utility_model = c("mixture", "wailoo"))
#' input.dat <- get_input_data(pop = pop)
#' parsamp <- sample_pars(n = 10, input_dat = input.dat)
#' sim.out <- sim_iviRA(tx_seqs = tx.seq, input_data = input.dat, pars = parsamp,
#'                     model_structures = mod.structs, output = "data")
#' head(sim.out)
#' 
#' @export 
sim_iviRA <- function(tx_seqs, input_data, pars, model_structures, 
                      max_months = NULL, tx_data = iviRA::treatments,
                      hist = c("naive", "experienced"),
                      output = c("data", "summary"), 
                      discount_qalys = .03, discount_cost = .03){
  hist <- match.arg(hist)
  output <- match.arg(output)
  
  # Prepping for the simulation
  ## check correct object types used as arguments
  if (!inherits(input_data, "input_data")){
    stop("The argument 'input_data' must be of class 'input_data'")
  }
  if (!inherits(model_structures, "model_structures")){
      stop("The argument 'model_structures' must be of class 'model_structures'")
  }
  if (!inherits(pars, "par_sample")){
    stop("The argument 'pars' must be of class 'par_sample'")
  }
  
  ## treatment indices
  if (is.vector(tx_seqs)) tx_seqs <- matrix(tx_seqs, nrow = 1)
  tx.inds <- matrix(match(tx_seqs, tx_data$sname), nrow = nrow(tx_seqs),
                    ncol = ncol(tx_seqs))
  if(any(is.na(tx.inds))){
    stop("At least one treatment in tx_seqs is not in tx_data.")
  }
  
  ## check for missing values
  if ("das28-switch" %in% model_structures[, "tx_iswitch"]){
    if (any(is.na(pars$das28$d[,,unique(tx.inds)]))){
      stop(paste0("Parameter values relating treatment directly to change in DAS28 at 6 months ",
                  "(i.e., from an NMA) are missing for ",
                  "at least one of the treatments in your treatment sequences."))
    }
  }
  
  if ("haq" %in% model_structures[, "tx_ihaq"]){
    if (any(is.na(pars$haq$d[,,unique(tx.inds)]))){
      stop(paste0("Parameter values relating treatment directly to change in HAQ at 6 months ",
            "(i.e., from an NMA) are missing for ",
            "at least one of the treatments in your treatment sequences."))
    }
  }
  
  ## indexing
  cdmards.ind <- which(iviRA::treatments$sname == "cdmards") - 1
  nbt.ind <- which(iviRA::treatments$sname == "nbt") - 1
  tx.inds <- tx.inds - 1
  
  ## default internal values
  treat_gap <- 0
  if (is.null(max_months)){
    max_months <- 12 * 150
  }
  prob.switch.da <- matrix(rep(c(0, 0, 1, 1), each = pars$n), ncol = 4)
  
  # convert factors variables to character variables in dataframes
  tx_data[, route := as.character(route)]

  ## survival parameters
  ttd.dist <- model_structures[1, "ttd_dist"]
  pars.ttd.all <- pars$ttd.all[[ttd.dist]]
  pars.ttd.da <- pars$ttd.da[[ttd.dist]]
  pars.ttd.em <- pars$ttd.eular$moderate[[ttd.dist]]
  pars.ttd.eg <- pars$ttd.eular$good[[ttd.dist]]
  si.dist <- "exp"
  pars.si <- pars$ttsi
  
  ## time to treatment discontinuation
  x.ttd.all <- input_data$x.ttd.all
  x.ttd.eular <- input_data$x.ttd.eular
  x.ttd.da <- input_data$x.ttd.da
  
  ## treatment costs
  tc <- pars$tx.cost
  tc.seqs <- cbind(tx_seqs, "nbt")
  lookup.inds <- match(tc.seqs, tc$lookup$sname)
  agents <- aperm(array(match(unlist(tc$lookup[lookup.inds, -1, with = FALSE]),
                         iviRA::tx.cost$cost$sname) - 1,
                   dim = c(nrow(tc.seqs), ncol(tc.seqs), ncol(tc$lookup) - 1)),
                  perm = c(2, 3, 1))
  
  ## utility parameters
  wailoo.coefnames <- c("int", "age", "dis_dur", "haq0", "male", "prev_dmards", "haq")
  wailoo.colindx <- match(wailoo.coefnames, colnames(pars$utility.wailoo))
  if(any(is.na(wailoo.colindx))) {
    stop("Matrix utility.wailoo in list 'pars' must have column names 'int', 'age',
         'dis_dur', 'haq0', 'male', 'prev_dmards', and 'haq'")
  }
  parsamp.utility.wailoo <- pars$utility.wailoo[, wailoo.colindx, drop = FALSE]

  # Running the simulation
  sim.out <- sim_iviRA_C(tx_inds = tx.inds, tx_data = tx_data,
                         model_structures_mat = model_structures, hist = hist,
                         haq0 = input_data$haq0, das28_0 = input_data$das28,
                     sdai0 = input_data$sdai, cdai0 = input_data$cdai,
                     age0 = input_data$age, male = input_data$male, 
                     prev_dmards = input_data$prev.dmards,
                     nma_acr_list = pars$acr,
                     x_acr = input_data$x.acr,
                     nma_haq_list = pars$haq,
                     x_haq = input_data$x.haq,
                     nma_das28_list = pars$das28,
                     x_das28 = input_data$x.das28,
                     acr2eular = pars$acr2eular, acr2haq = pars$acr2haq, eular2haq = pars$eular2haq,
                     acr2das28 = pars$acr2das28, acr2sdai = pars$acr2sdai, acr2cdai = pars$acr2cdai,
                     tswitch_da = prob.switch.da,
                     haq_lprog_therapy = pars$haq.lprog.tx, haq_lprog_age = pars$haq.lprog.age,
                     haq_lcgm_delta = pars$haq.lcgm$delta, haq_lcgm_beta = pars$haq.lcgm$beta, 
                     rebound_factor = pars$rebound, 
                     lifetable_male = pars$lt$male, lifetable_female = pars$lt$female,
                     x_mort = input_data$x.mort, logor_mort = pars$mort.logor, 
                     x_ttd_all = x.ttd.all, x_ttd_da = x.ttd.da, x_ttd_eular = x.ttd.eular,
                     ttd_all_list = pars$ttd.all, ttd_da_list = pars$ttd.da,
                     ttd_eular_mod_list = pars$ttd.eular$moderate, ttd_eular_good_list = pars$ttd.eular$good,
                     cdmards = cdmards.ind, nbt = nbt.ind,
                     si_loc = pars.si, 
                     si_anc1 = matrix(NA, nrow = pars$n, ncol = ncol(pars.si)), 
                     si_anc2 = matrix(NA, nrow = pars$n, ncol = ncol(pars.si)), 
                     si_dist = si.dist, 
                     haqdelta_loghr = pars$mort.loghr.haqdif, max_months = max_months,
                     hosp_days = pars$hosp.cost$hosp.days, cost_pday = pars$hosp.cost$cost.pday,
                     mgmt_cost = rowSums(pars$mgmt.cost), 
                     si_cost = pars$si.cost, prod_loss = pars$prod.loss,
                     tc_list = c(list(agents = agents), tc[c("cost", "discount")]), 
                     weight = input_data$weight, 
                     coefs_wailoo = parsamp.utility.wailoo, 
                     pars_util_mix = pars$utility.mixture, si_ul = pars$si.ul,
                     utility_tx_attr = list(data = input_data$x.attr, 
                                    pars = pars$utility.tx.attr),
                     discount_rate = list(qalys = discount_qalys, cost = discount_cost),
                     output = output)
  if (output == "data"){
      sim.out <- as.data.table(sim.out)
      
      ## C++ to R indices
      sim.out[, model := model + 1]
      sim.out[, sim := sim + 1]
      sim.out[, id := id + 1]
      sim.out[, line := line + 1]
      sim.out[, tx_cycle := tx_cycle + 1]
      sim.out[, tx := tx + 1]
  } else{
      names(sim.out)[names(sim.out) == "means1"] <- "means"
      sim.out$means <- data.table(cbind(sim.out$means, sim.out$means2))
      sim.out$means2 <- NULL
      sim.out$time.means <- data.table(sim.out$time.means)
      sim.out$time.means <- sim.out$time.means[alive > 0]
      sim.out$out0 <- data.table(sim.out$out0)
    
      ## C++ to R indices
      sim.out$means[, model := model + 1]
      sim.out$means[, sim := sim + 1]
      sim.out$time.means[, model := model + 1]
      sim.out$time.means[, sim := sim + 1]
      sim.out$out0[, model := model + 1]
      sim.out$out0[, sim := sim + 1]
      sim.out$out0[, id := id + 1]
      sim.out$out0[, tx := tx + 1]
  }
  
  # Return
  return(sim.out)
}

#' Simulate utility after IPS
#'
#' Simulate utility after running \link{sim_iviRA} with \code{output = "data"}. This can be useful 
#' in cases where you want to use a different algorithm to estimate utility, but do not want to rerun 
#' the entire simulation.
#' @param simhaq Simulation output from \link{sim_iviRA}. Must include columns \code{yrlen} for
#' year length of model cycle, \code{sim} for simulation number, and \code{si} for whether a serious
#' infection occured during the model cycle. 
#' @param male Indicator = 1 for males and 0 for females.
#' 
#' @details Note that disease duration is set to 18.65 years in \code{sim_utility_wailoo}, which
#' is the mean value from the Wailoo (2006) paper used for the parameter estimates. Age and 
#' the HAQ score are taken from the simulation output.
#' 
#' @return For \code{sim_utility_mixture} and \code{sim_utility_wailoo}, a vector of 
#' simulated utility for each row returned in \code{simhaq}. For \code{sim_qalys}, a vector
#' of QALYs for each row in \code{simhaq}.
#' 
#' @examples
#' pop <- sample_pop(n = 10)
#' tx.seq <- c("adamtx", "cdmards")
#' mod.structs <- select_model_structures(utility_model = "wailoo")
#' input.dat <- get_input_data(pop = pop)
#' parsamp <- sample_pars(n = 10, input_dat = input.dat)
#' sim.out <- sim_iviRA(tx_seqs = tx.seq, input_data = input.dat, pars = parsamp,
#'                     model_structures = mod.structs, output = "data")
#' utility.mix <- sim_utility_mixture(simhaq = sim.out, male = pop[, "male"], pars = parsamp$utility.mixture)
#' utility.wailoo <- sim_utility_wailoo(simhaq = sim.out, haq0 = pop[, "haq0"], male = pop[, "male"],
#'                                 prev_dmards = pop[, "prev_dmards"], 
#'                                 coefs = parsamp$utility.wailoo)
#' qalys.mix <- sim_qalys(simhaq = sim.out, utility = utility.mix, si_ul = parsamp$si.ul,
#'                        x_attr = input.dat$x.attr, tx_attr_coef = parsamp$utility.tx.attr)
#' head(utility.mix)
#' head(utility.wailoo)   
#' head(qalys.mix)
#' 
#' @name sim_utility
NULL
#> NULL

#' @param pars List of sampled parameters needed to simulate utility using the Hernandez Alava (2013) mixture 
#' model (i.e., the element \code{utility.mixture} returned by \link{sample_pars}).
#' to \link{sample_pars}.
#' 
#' @rdname sim_utility
#' @export
sim_utility_mixture <- function(simhaq, male, pars){
  check_sim_utility_mixture(simhaq, male, pars)
  util <- iviRA:::sim_utility_mixtureC(id = simhaq$id - 1, sim = simhaq$sim - 1, haq = simhaq$haq, 
                               pain_mean = pars$pain$pain.mean, haq_mean = pars$pain$haq.mean, 
                               pain_var = pars$pain$pain.var, haq_var = pars$pain$haq.var, 
                               painhaq_cor = pars$pain$painhaq.cor, 
                               age = simhaq$age, male = male,
                               beta1 = pars$beta1, beta2 = pars$beta2, beta3 = pars$beta3,
                               beta4 = pars$beta4, 
                               alpha1 = pars$alpha1, alpha2 = pars$alpha2, alpha3 = pars$alpha3, 
                               alpha4 = pars$alpha4,
                               alpha = pars$alpha,
                               epsilon1 = pars$epsilon1, epsilon2 = pars$epsilon2, epsilon3 = pars$epsilon3, 
                               epsilon4 = pars$epsilon4,
                               mu = pars$mu, delta = pars$delta)
  return(util)
}

#' Check parameters of sim_utility_mixture
#'
#' Error messages when incorrect inputs are passed to sim_utility_mixture.
#' 
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

#' @param simhaq Simulation output from \link{sim_iviRA}. Variables needed are \code{sim}, \code{id},
#'  \code{age}, and \code{haq}.
#' @param haq0 HAQ score at baseline.
#' @param prev_dmards Number of previous DMARDs.
#' @param coefs Matrix of sampled coefficients needed to simulate utility using the Wailoo (2006) model 
#' (i.e., the element \code{utility.wailoo} returned by \link{sample_pars}). Note that the matrix 
#' columns must contain the exact same variables as generated by \link{sample_pars}. 
#' 
#' @rdname sim_utility
#' @export
sim_utility_wailoo <- function(simhaq, haq0, male, prev_dmards,
                               coefs){
  check_sim_utility_wailoo(simhaq, haq0, male, prev_dmards, coefs)
  dis.dur <- 18.65 # based on mean disease duration in Wailoo (2006)
  util <- iviRA:::sim_utility_wailooC(sim = simhaq$sim - 1, id = simhaq$id - 1, age = simhaq$age,
                              disease_duration = dis.dur, haq0 = haq0, male = male, 
                              prev_dmards = prev_dmards, haq = simhaq$haq,
                              b_int = coefs[, "int"], b_age = coefs[, "age"], b_disease_duration = coefs[, "dis_dur"], 
                              b_haq0 = coefs[, "haq0"], b_male = coefs[, "male"],
                              b_prev_dmards = coefs[, "prev_dmards"], b_haq = coefs[, "haq"])
  return(util)
}

#' Check parameters of sim_utility_wailoo
#'
#' Error messages when incorrect inputs are passed to sim_utility_wailoo.
#' 
#' @param simhaq Simulation output from \link{sim_iviRA}
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
  if(ncol(coefs) != 7) stop("Number of columns in 'coefs' must be 7")
  if(nrow(coefs) != n) stop(paste0("Number of rows in 'coefs' must be equal to the number of ",
                           "sampled parameter sets which equals ", n))
  if(any(!colnames(coefs) %in% c("int", "age", "dis_dur", "haq0", "male", "prev_dmards", "haq"))){
    stop(paste0("Column names of 'coefs' must be (in order) 'int', 'age', 'dis_dur'",
                "'haq0', 'male', 'prev_dmards', and 'haq'."))
  }
}

#' @param utility Simulated utility from \code{sim_utility_mixture} or \code{sim_utility_wailoo}.
#' @param si_ul Sampled utility loss. Equivalent to output \code{si.ul} from \link{sample_pars}.
#' @param x_attr Treatment attribute data (e.g., ouptut \code{x.attr} from
#' \link{get_input_data})
#' @param tx_attr_coef Distribution of coefficient vector for treatment attributes (e.g., output
#' \code{utility.tx.attr} from \link{sample_pars}.)
#' 
#' @rdname sim_utility
#' @export
sim_qalys <- function(simhaq, utility, si_ul, x_attr, tx_attr_coef){
  check_sim_qalys(simhaq, utility, si_ul, x_attr, tx_attr_coef)
  qalys <- sim_qalysC(utility = utility, yrlen = simhaq$yrlen, 
                      sim = simhaq$sim - 1, tx = simhaq$tx - 1,
                      si = simhaq$si, si_ul = si_ul,
                      x_attr = x_attr, tx_attr_coef = tx_attr_coef)
  return(qalys)
}

#' Check parameters of sim_qalys
#'
#' Error messages when incorrect inputs are passed to sim_qalys
#' 
#' @param si_ul \code{si_ul} as passed to \link{sim_qalys}
#' @keywords internal
check_sim_qalys <- function(simhaq, utility, si_ul, x_attr, tx_attr_ug){
  # utility 
  if(is.null(utility)) stop("'utility' is missing")
  
  # simhaq
  if(is.null(simhaq$yrlen)) stop("'yrlen' column of 'simhaq' is missing")
  if(is.null(simhaq$sim)) stop("'sim' column of 'simhaq' is missing")
  if(is.null(simhaq$si)) stop("'si' column of 'simhaq' is missing")
  if(is.null(simhaq$tx)) stop("'tx' column of 'simhaq' is missing")
  
  # si_ul
  n <- max(simhaq$sim)
  if(is.null(si_ul)) stop("'si_ul' is missing")
  if(!is.vector(si_ul)) stop("'si_ul' must be a vector")
  if(length(si_ul) != n){
    stop(paste0("Number of rows in 'si_ul' must ",
                "be equal to the number of sampled parameter sets which equals ", n))
  }
  
  # x_attr
  if(is.null(x_attr)) stop("'x_attr' is missing")
  if(!is.matrix(x_attr)) stop("'x_attr' must be a matrix")
  
  # tx_attr_ug
  if(is.null(tx_attr_ug)) stop("'tx_attr_ug' is missing")
  if(!is.matrix(tx_attr_ug)) stop("'tx_attr_ug' must be a matrix")
}

#' Simulate change in HAQ at 6 months
#'
#' Simulate change in HAQ at 6 months using different model structures relating treatment to HAQ.
#' 
#' @param treatments Vector of treatments to simulate. Should correpond to \code{sname} in 
#' iviRA::treatments.
#' @param input_data An object of class \code{input_data} returned from \link{get_input_data}. 
#' The only elements required are \code{x_acr} and \code{x_haq}.
#' @param pars List of sampled parameter values generated from \link{get_input_data}. The only 
#' elements required are \code{acr}, \code{haq}, \code{acr2haq}, \code{acr2eular}, and
#' \code{eular2haq}.
#' @param tx_lookup Vector of treatments with names equivalent to \code{iviRA::treatments$sname}.
#' Index of treatments in \code{treatments} are matched against treatments in \code{tx_lookup} by 
#' name. Indices of parameter estimates from the network meta-analyses bust be in the same order as
#' in \code{tx_data}.
#' @param hist Is the patient tDMARD naive or tDMARD experienced? 
#' @param line Line of therapy
#' @param tx_ihaq Equivalent to argument \code{tx_haq} in \link{select_model_structures}.
#' 
#' @return A \code{data.table} with the following columns: 
#' \describe{
#' \item{tx}{Treatment index based on \code{tx_lookup}.}
#' \item{model}{The model structure component relating treatment to change in HAQ 
#' (i.e., \code{tx_ihaq}).}
#' \item{sim}{Simulation number denoting a randomly sampled parameter set from \link{sample_pars}.}
#' \item{id}{ID number denoting a simulated patients (e.g., from \link{sample_pop}).}
#' \item{dhaq}{Change in HAQ score from baseline at 6 months.}
#' }
#'
#' @examples 
#' pop <- sample_pop(n = 10)
#' input.dat <- get_input_data(pop = pop)
#' parsamp <- sample_pars(n = 10, input_data = input.dat)
#' sim <- sim_dhaq6(treatments = c("cdmards", "adamtx"), input_data = input.dat,
#'                  pars = parsamp, tx_ihaq = c("acr-haq", "acr-eular-haq"))
#' head(sim)                 
#' tail(sim)
#' 
#' @export
sim_dhaq6 <- function(treatments, input_data, pars,
                      tx_lookup = iviRA::treatments$sname,
                      hist = c("naive", "experienced"), 
                      line = 1, 
                      tx_ihaq = c("acr-haq", "acr-eular-haq", "haq")){
  hist <- match.arg(hist)
  nbt.ind <- which(tx_lookup == "nbt") - 1
  cdmards.ind <- which(tx_lookup == "cdmards") - 1
  treatments.ind <- match(treatments, tx_lookup) - 1
  sim <- sim_dhaq6C(npats = input_data$n, nsims = pars$n, 
                    hist = hist, line = line - 1,
                    tx_inds = treatments.ind, nbt_ind = nbt.ind,
                    x_acr = input_data$x.acr, nma_acr_list = pars$acr,
                    x_haq = input_data$x.haq, nma_haq_list = pars$haq,
                    acr2eular = pars$acr2eular, acr2haq = pars$acr2haq,
                    eular2haq = pars$eular2haq, tx_ihaq_type = tx_ihaq)
  sim <- data.table(sim)
  sim[, tx := tx + 1]
  sim[, sim := sim + 1]
  sim[, id := id + 1]
  return(sim)
}