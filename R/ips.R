#' The IVI-RA simulation
#'
#' Run the IVI-RA individual patient simulation model.
#' 
#' @param arms Name of arms in the treatment treatment sequence. May be a vector consisting of 
#' a single treatment sequence or a matrix of unique sequences for each patient.
#' @param input_data An object of class 'input_data' returned from \link{get_input_data}.
#' @param pars List of sampled parameter values generated from \link{sample_pars}.
#' @param model_structure Object of class model structure generated from \link{select_model_structure}.
#' @param max_months Maximum number of months to run the model for. Default is NULL which implies that
#' the model is simulated over each patient's lifetime.
#' @param treatments_lookup Vector of names of all treatments included in the parameter
#' estimates. Index of of treatments in \code{arms} are matched against treatments in
#' \code{treatments_lookup} by name. Indices of treatment-specific parameter estimates must be 
#' in the same order as treatments in \code{treatments_lookup}.   
#' @param output Specifies the format of output returned from the simulation. Options are \code{data} 
#' and \code{summary}. When \code{data} is specified, each simulated value (i.e, by model,
#' sampled parameter set, individual, and time-period) is returned in a \code{data.table}. If 
#' \code{summary} is selected, then only summary measures are returned.
#' @param discount_qalys Discount rate for QALYs. Only used when \code{output = "summary"}; 
#' otherwise, discounts can be applied to the simulated output.
#' @param discount_cost Discount rate for cost variables. Only used when \code{output = "summary"};
#' otherwise, discounts can be applied to the simulated output.
#' 
#' @return 
#' The \code{output = "data"} options returns all simulated output. However, since output is 
#' returned for each model, sampled parameter set, individual, and model cycle, 
#' the size of the output can be very large and can cause memory management issues. The 
#' \code{output = "summary"}, which only provides summaries of the simulation output, can be 
#' useful in the cases. 
#' 
#' When \code{output = "data"} is selected, the simulation returns a \code{data.table} with the 
#' following columns:
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
#' \item{si_cost}{Cost due to serious infections.}
#' \item{prod_loss}{Productivity loss (i.e., lost earnings).}
#' \item{utility}{Simulated utility score.}
#' \item{qalys}{Quality-adjusted life-years (QALYs).}
#'}
#'
#' @export 
sim_iviRA <- function(arms, input_data, pars, model_structures, 
                      max_months = NULL, treatment_lookup = iviRA::treatments$sname,
                      output = c("data", "summary"), 
                      discount_qalys = .03, discount_cost = .03,
                      check = TRUE){
  output <- match.arg(output)
  
  # PREPPING FOR THE SIMULATION
  ## check correct object types used as arguments
  if (!inherits(input_data, "input_data")){
    stop("The argument 'input_data' must be of class 'input_data'")
  }
  if (!inherits(model_structures, "model_structures")){
      stop("The argument 'model_structures' must be of class 'model_structures'")
  }
  
  ## treatment arm indices
  if (is.vector(arms)) arms <- matrix(arms, nrow = 1)
  arminds <- matrix(match(arms, treatment_lookup), nrow = nrow(arms),
                    ncol = ncol(arms))
  
  ## indexing
  cdmards.ind <- which(iviRA::treatments$sname == "cdmards") - 1
  nbt.ind <- which(iviRA::treatments$sname == "nbt") - 1
  arminds <- arminds - 1
  
  ## default internal values
  treat_gap <- 0
  if (is.null(max_months)){
    max_months <- 12 * 150
  }
  prob.switch.da <- matrix(rep(c(0, 0, 1, 1), each = pars$n), ncol = 4)

  ## survival parameters
  ttd.dist <- model_structures[1, "ttd_dist"]
  pars.ttd.all <- pars$ttd.all[[ttd.dist]]
  pars.ttd.da <- pars$ttd.da[[ttd.dist]]
  pars.ttd.em <- pars$ttd.eular$moderate[[ttd.dist]]
  pars.ttd.eg <- pars$ttd.eular$good[[ttd.dist]]
  si.dist <- "exp"
  pars.si <- pars$ttsi
  
  ## time to treatment discontinuation
  if (is.null(input_data$x.ttd.all)){
    x.ttd.all <- matrix()
  } else{
    x.ttd.all <- input_data$x.ttd.all
  }
  
  if (is.null(input_data$x.ttd.eular)){
    x.ttd.eular <- matrix()
  } else{
    x.ttd.eular <- input_data$x.ttd.eular
  }
  
  if (is.null(input_data$x.ttd.da)){
    x.ttd.da <- matrix()
  } else{
    x.ttd.da <- input_data$x.ttd.da
  }
  
  ## treatment costs
  tc <- pars$treat.cost
  tc.arms <- cbind(arms, "nbt")
  lookup.inds <- match(tc.arms, tc$lookup$sname)
  agents <- aperm(array(match(unlist(tc$lookup[lookup.inds, -1, with = FALSE]),
                         iviRA::treat.cost$cost$sname) - 1,
                   dim = c(nrow(tc.arms), ncol(tc.arms), ncol(tc$lookup) - 1)),
                  perm = c(2, 3, 1))
  
  ## utility parameters
  wailoo.coefnames <- c("int", "age", "dis_dur", "haq0", "male", "prev_dmards", "haq")
  wailoo.colindx <- match(wailoo.coefnames, colnames(pars$utility.wailoo))
  if(any(is.na(wailoo.colindx))) {
    stop("Matrix utility.wailoo in list 'pars' must have column names 'int', 'age',
         'dis_dur', 'haq0', 'male', 'prev_dmards', and 'haq'")
  }
  parsamp.utility.wailoo <- pars$utility.wailoo[, wailoo.colindx, drop = FALSE]

  # RUNNING THE SIMULATION
  sim.out <- sim_iviRA_C(arm_inds = arminds, model_structures_mat = model_structures,
                         haq0 = input_data$haq0, das28_0 = input_data$das28,
                     sdai0 = input_data$sdai, cdai0 = input_data$cdai,
                     age0 = input_data$age, male = input_data$male, 
                     prev_dmards = input_data$prev.dmards,
                     nma_acr1 = pars$acr$p1, nma_acr2 = pars$acr$p2, 
                     nma_dhaq1 = pars$haq$dy1, nma_dhaq2 = pars$haq$dy2,
                     nma_das28_1 = pars$das28$dy1, nma_das28_2 = pars$das28$dy2,
                     acr2eular = pars$acr2eular, acr2haq = pars$acr2haq, eular2haq = pars$eular2haq,
                     acr2das28 = pars$acr2das28, acr2sdai = pars$acr2sdai, acr2cdai = pars$acr2cdai,
                     tswitch_da = prob.switch.da,
                     haq_lprog_therapy = pars$haq.lprog.tx, haq_lprog_age = pars$haq.lprog.age,
                     haq_lcgm_delta = pars$haq.lcgm$delta, haq_lcgm_beta = pars$haq.lcgm$beta, 
                     rebound_factor = pars$rebound, 
                     lifetable_male = pars$lt$male, lifetable_female = pars$lt$female,
                     x_mort = input_data$x.mort, logor_mort = pars$logor.mort, 
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
                     discount_rate = list(qalys = discount_qalys, cost = discount_cost),
                     output = output)
  if (output == "data"){
      sim.out <- as.data.table(sim.out)
      
      ## C++ to R indices
      sim.out[, model := model + 1]
      sim.out[, sim := sim + 1]
      sim.out[, id := id + 1]
      sim.out[, tx_seq := tx_seq + 1]
      sim.out[, tx_cycle := tx_cycle + 1]
      sim.out[, tx := tx + 1]
  } else{
      sim.out$means <- data.table(sim.out$means)
      sim.out$time.means <- data.table(sim.out$time.means)
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
  
  # RETURN
  return(sim.out)
}

#' Select model structures
#'
#' Select the model structures to be used in the IVI-RA individual patient simulation.
#' 
#' @param tx_ihaq Model structure relating treatment to HAQ during the first 6 months of 
#' treatment. Options, which are equivalent to H1-H3 in the documentation are:
#' \itemize{
#' \item{\code{acr-haq}}{ H1: Treatment -> ACR -> HAQ }
#' \item{\code{acr-eular-haq}}{ H2: Treatment -> ACR -> EULAR -> HAQ}
#' \item{\code{haq}}{ H3:Treatment -> HAQ}
#' }
#' @param tx_iswitch Model structure relating treatment to switching during the first
#' 6 months of treatment. Options, which are equivalent to S1-S6 in the documentation are:
#' \itemize{
#' \item{\code{acr-switch}}{ S1: Treatment -> ACR -> Switch}
#' \item{\code{acr-das28-switch}}{ S2: Treatment -> ACR -> DAS28 -> Switch}
#' \item{\code{acr-sdai-switch}}{ S3: Treatment -> ACR -> SDAI -> Switch}
#' \item{\code{acr-cdai-switch}}{ S4: Treatment -> ACR -> CDAI -> Switch}
#' \item{\code{das28-switch}}{ S5: Treatment -> DAS28 -> Switch}
#'  \item{\code{acr-eular-switch}}{ S6: Treatment -> ACR -> EULAR -> Switch}
#' 
#' }
#' @param cdmards_haq_model Model used for long-term HAQ progression. Options are:
#' \itemize{
#' \item{lcgm}{ Latent class growth model}
#' \item{linear}{ Constant linear HAQ progression}
#' }
#' If \code{lgcm} is chosen, then a latent class growth model is used for cDMARDs
#' and NBT but a constant annual rate is is assumed for all other therapies; otherwise 
#' a constant linear HAQ progression is assumed for all therapies including cDMARDs and NBT.
#' @param ttd_dist Distribution used to model time to treatment discontinuaton. Options are:
#' \itemize{
#' \item{exponential}{ Exponential}
#' \item{weibull}{ Weibull}
#' \item{gompertz}{ Gompertz}
#' \item{gamma}{ Gamma}
#' \item{llogis}{ Log-logistic}
#' \item{lnorm}{ Lognormal}
#' \item{gengamma}{ Generalized gamma}
#' }
#' @param utility_model Model used to estimate patient utility as a function of HAQ and patient
#' characteristics. Options are:
#' \itemize{
#' \item{mixture}{ Hernandez Alava (2013) mixutre model}
#' \item{wailoo}{ Wailoo 2006 logistc regression}
#' }
#' @export
select_model_structures <- function(tx_ihaq = "acr-haq",
                                   tx_iswitch = "acr-switch",
                                   cdmards_haq_model = "lcgm", 
                                   ttd_dist = "exponential",
                                   utility_model = "mixture"){
  # 
  n <- vector(length = 5)
  n[1] <- length(tx_ihaq)
  n[2] <- length(tx_iswitch)
  n[3] <- length(cdmards_haq_model)
  n[4] <- length(ttd_dist)
  n[5] <- length(utility_model)
  if (max(n) > 1){
    n.g1 <- n[n > 1]
    max.n.g1 <- max(n.g1)
    if (any(max.n.g1 - n.g1 > 0)){
      stop("Length of all vectors must be the same")
    }
    if (any(max.n.g1 - n > 0)){
      if (n[1] == 1){
        tx_ihaq <- rep(tx_ihaq, max.n.g1)
      }
      if (n[2] == 1){
        tx_iswitch <- rep(tx_iswitch, max.n.g1)
      }
      if (n[3] == 1){
        cdmards_haq_model <- rep(cdmards_haq_model, max.n.g1)
      }
      if (n[4] ==1){
        ttd_dist <- rep(ttd_dist, max.n.g1)
      } 
      if (n[5] ==1){
        utility_model <- rep(utility_model, max.n.g1)
      } 
    }
  }
  
  # are valid options selected?
  if (any(!tx_ihaq %in% c("acr-haq", "acr-eular-haq", "haq"))){
      stop("Values in 'tx_ihaq' must be 'acr-haq', 'acr-eular-haq' or 'haq'.")
  }
  
  if (any(!tx_iswitch %in% c("acr-switch", "acr-das28-switch",
                             "acr-sdai-switch", "acr-cdai-switch", 
                             "das28-switch", "acr-eular-switch"))){
      stop(paste0("Values in 'tx_iswitch' must be 'acr-switch', 'acr-das28-switch', 'acr-sdai-switch',",
                  " 'cr-cdai-switch', 'das28-switch', or 'acr-eular-switch'."))
  }
  
  if (any(!cdmards_haq_model %in% c("lcgm", "linear"))){
      stop("Values in 'cdmards_haq_model' must be 'lcgm' or 'linear'.")
  } 
  
  if (any(!ttd_dist %in% c("exponential", "weibull", "gompertz", 
                           "gamma", "llogis", "lnorm", "gengamma"))){
    stop(paste0("Values in 'ttd_dist' must be 'exponential', 'weibull', 'gompertz', 'gamma',", 
         " 'llogis', 'lnorm', or 'gengamma'."))
  } 

  if (any(!utility_model %in% c("mixture", "wailoo"))){
      stop("Values in 'utility_model' must be 'mixture' or 'wailoo'.")
  } 
  
  # are valid combinations of options selected?
  ## tx_ihaq = acr-haq
  val <- ifelse(tx_ihaq == "acr-haq" & tx_iswitch == "acr-eular-switch", 1, 0)
  if (any(val > 0)){
    stop("'tx_iswitch' option 'acr-eular-switch' cannot be used with 'tx_ihaq' option
         'acr-haq'.")
  }
  
  ## tx_ihaq = haq
  val <- ifelse(tx_ihaq == "haq" & tx_iswitch != "das28-switch", 1, 0)
  if (any(val == 1)){
    stop("When 'tx_ihaq' option 'haq' is selected, 'tx_iswitch' must equal 'das28-switch'.")
  }

  # return
  model.structure <- matrix(c(tx_ihaq, tx_iswitch, cdmards_haq_model, ttd_dist, utility_model), ncol = 5)
  colnames(model.structure) <- c("tx_ihaq", "tx_iswitch", "cdmards_haq_model", "ttd_dist", "utility_model")
  class(model.structure) <- "model_structures"
  return(model.structure)
}

#' Check parameters for ivi_RA
#'
#' Error messages when incorrect parameters are passed to ivi_RA.
#' 
#' @param input_data \code{input_data} as passed to \link{sim_haq}.
#' @param pars \code{pars} as passed to \link{sim_haq}.
#' @param tx_ihaq Treatment to HAQ pathway first 6 months.
#' @param tx_iswitch Treatment switching pathway first 6 months.
#' @keywords internal
check_pars <- function(arminds, input_data, pars, tx_ihaq, tx_iswitch){
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
  if (tx_ihaq %in% c("acr-haq", "acr-eular-haq")){
    if (any(is.na(nma.acr1) | is.na(nma.acr2))){
      stop ("'acr' element of pars has missing parameter values for one of the selected treatment 
            arms; that is, the NMA results are missing.")
    }
  }

  ## NMA HAQ
  nma.haq1 <- pars$haq$dy1[, arminds.unique]
  nma.haq2 <- pars$haq$dy2[, arminds.unique] 
  if (tx_ihaq == "haq"){
    if (any(is.na(nma.haq1) | is.na(nma.haq2))){
      stop ("'haq' element of pars has missing parameter values for one of the selected treatment 
            arms; that is, the NMA results are missing.")
    }
  }
  
  ## NMA DA28
  nma.das28.1 <- pars$das28$dy1[, arminds.unique]
  nma.das28.2 <- pars$das28$dy2[, arminds.unique] 
  if (tx_iswitch == "das28-switch"){
    if (any(is.na(nma.das28.1) | is.na(nma.das28.1))){
      stop ("'das28' element of pars has missing parameter values for one of the selected treatment 
            arms; that is, the NMA results are missing.")
    }
  }
}


#' Simulate utility using Hernandez-Alava mixture model
#'
#' Simulate utility from HAQ score simulated using using mixture model from
#' Hernandez-Alva (2013). Can be used to simulate utility after running \code{ivi_RA}.
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
  util <- sim_utility_mixtureC(id = simhaq$id - 1, sim = simhaq$sim - 1, haq = simhaq$haq, 
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

#' Simulate utility using Wailoo (2006) model
#'
#' Simulate utility from simulated HAQ score simulated using using model from
#' Wailoo (2006). Can be used to simulate utility after running \code{sim_iviRA}.
#' 
#' @param simhaq Simulation output from \link{sim_iviRA}. Variables needed are \code{sim}, \code{id},
#'  \code{age}, and \code{haq}.
#' @param haq0 HAQ score at baseline.
#' @param male Indicator = 1 for males and 0 for females.
#' @param prev_dmards Number of previous DMARDs.
#' @param coefs Matrix of coefficients needed to simulate utility using the Wailoo (2006) model. 
#' See the documentation in \link{sample_pars} for details. Note that the matrix columns must contain the
#' same variables as generated by \link{sample_pars}. 
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
  dis.dur <- 18.65 # based on mean disease duration in Wailoo (2006)
  util <- sim_utility_wailooC(sim = simhaq$sim - 1, id = simhaq$id - 1, age = simhaq$age,
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


