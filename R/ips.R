#' Simulate HAQ over time
#'
#' An individual patient simulation of HAQ scores for rheumatoid arthritis patients given a sequence of J 
#' treatments. An "outer" loop iterates over S draws from sampled parameter sets and an "inner" loop iterates of
#' N patients. Cycle lengths are 6 months long. The simulation is written in C++ for speed. 
#' 
#' @param arminds Indices of treatment arms consisting of sequences of therapies (each element
#' is the index of a therapy in \code{therapy.pars$info}. May be a vector consisting of a single treatment sequence or a matrix 
#' of unique sequences for each patient.
#' @param input_data List of input data. Required inputs are \code{haq0}, \code{age}, \code{male}, \code{x.mort}, 
#' and \code{x.dur} as generated from \link{input_data}.
#' @param pars List of parameters. Required parameters are \code{rebound}, \code{acr1}, \code{acr2}, \code{acr2eular}, \code{eular2haq}, 
#' \code{haq.lprog.therapy}, \code{haq.lprog.age}, \code{logor.mort}, \code{mort.loghr.haqdif}, \code{si.surv},
#' \code{ttd.eular.mod}, \code{ttd.eular.good
#' }, and \code{lt} as
#' generated from \link{sample_pars}.
#' @param dur_dist Survival distribution for treatment duration.
#' @param si_dist Survival distribution for serious infections. Currently rate is assumed
#' to be constant so an exponential distribution is used.
#' @param max_months Maximum number of months to run the model for. Default is NULL which implies that
#' the model is simulated over each patient's lifetime.
#' @param nbt_ind Index for the non-biologic (i.e. a term used to define a selection of treatments clinicians use
#' after the last biologic in a treatment sequence).
#' @param check Should the function check parameters and input data passed to the model? Default is TRUE.
#' 
#' @return A \code{data.table} with the following columns:
#'\describe{
#' \item{sim}{Simulation number. Indexed from 1 to S where S is the number of randomly sampled parameter sets (e.g. n from \link{sample_pars}).}
#' \item{id}{ID number. Indexed from 1 to N where N is the number of simulated patients (e.g. n from \link{sample_pats}).}
#' \item{month}{Month since a simulated patient began the first therapy in a treatment sequence.}
#' \item{therapy}{Index of therapy used. Given J total therapies, the first J - 1 therapies match the indices from \code{arminds}. The final \code{therapy}
#'  is always equal to the index of the non-biologic therapy (\code{nbt_ind}).}
#' \item{therapy_seq}{Number of therapies in a sequence. First therapy equal to 1, second therapy equal to 2, \ldots}
#' \item{therapy_cycle}{Number of model cycles since a patient began taking a given therapy in a treatment sequence. \code{therapy_cycle} = 1
#' during the initial 6-month treatment period.}
#' \item{death}{Equal to 1 if patient died during the model cycle and 0 otherwise. }
#' \item{age}{Age of patient which increases with the model cycles.}
#' \item{ttd}{Time to treatment discontinuation. Measured in terms of model cycles (e.g. ttd = 2 if therapy will discontinue in 
#' 1 year given 6-months cycles). \code{ttd} is measured at the end of each cycle. Patients switch therapies during the cycle in which \code{ttd} 
#' becomes negative and HAQ rebounds during the cycle.}
#' \item{acr}{Simulated ACR response during the initial 6-month period for a new therapy. Constant within \code{therapy}. 
#' Categories are 0 (ACR < 20), 1 (ACR 20-50), 2 (ACR 50-70), and 3 (ACR 70+).}
#' \item{eular}{Simulated EULAR response during the initial 6-month period for a new therapy. Constant within \code{therapy}. 
#' Categories are 0 (no EULAR response), 1 (moderate EULAR response), and 2 (good EULAR response).}
#' \item{haq}{HAQ score. Restricted to range between 0 and 3.}
#' \item{ttsi}{Time to serious infection. Like \code{ttd}, measured in terms of model cycles. \code{ttsi} is measured at the end of each 
#' cycle.}
#' \item{si}{Equal to 1 if treatment discontinuation was caused by a serious infection and 0 otherwise.}
#' \item{yrlen}{Length of a model cycle in years. Equal to 0.5 given 6-month cycles.}
#'}
#' @export
sim_haq <- function(arminds, input_data, pars, dur_dist = "lnorm", si_dist = "exp",
                   max_months = NULL, nbt_ind = which(therapy.pars$info$sname == "nbt"), check = TRUE){
  if (check) check_sim_haq(input_data, pars)
  treat_gap <- 0
  cycle_length <- 6
  nbt.ind <- nbt_ind - 1
  arminds <- arminds - 1
  if (class(arminds) == "numeric") arminds <- matrix(arminds, nrow = 1)
  if (!nrow(arminds) %in% c(1, length(input_data$haq0))){
    stop("Number of treatment sequences must either be the same for each patient or equal to the number of patients
         in input_data")
  }
  if (is.null(max_months)){
    max_months <- 12 * 150
  }
  pars.ttd.em <- pars$ttd.eular.mod[[dur_dist]]
  pars.ttd.eg <- pars$ttd.eular.good[[dur_dist]] 
  pars.si <- pars$si.surv[[si_dist]]
  simout <- sim_haqC(arminds, input_data$haq0, input_data$age, 
                  input_data$male,
                  pars$acr$p1, pars$acr$p2, pars$acr2eular, pars$eular2haq, 
                  pars$haq.lprog.therapy, pars$haq.lprog.age,
                  pars$rebound, pars$lt$male, pars$lt$female,
                 input_data$x.mort, pars$logor, dur_dist, input_data$x.dur, 
                 pars.ttd.em$sample[, pars.ttd.em$loc.index, drop = FALSE], 
                 pars.ttd.em$sample[, pars.ttd.em$anc1.index, drop = FALSE],
                 pars.ttd.em$sample[, pars.ttd.em$anc2.index, drop = FALSE], 
                 pars.ttd.eg$sample[, pars.ttd.eg$loc.index, drop = FALSE], 
                 pars.ttd.eg$sample[, pars.ttd.eg$anc1.index, drop = FALSE],
                 pars.ttd.eg$sample[, pars.ttd.eg$anc2.index, drop = FALSE],
                 cycle_length, treat_gap, nbt.ind,
                 pars.si$sample[, pars.si$loc.index, drop = FALSE], 
                 pars.si$sample[, pars.si$anc1.index, drop = FALSE],
                 pars.si$sample[, pars.si$anc2.index, drop = FALSE], 
                 si_dist, pars$mort.loghr.haqdif, max_months)
  simout <- as.data.table(simout)
  colnames(simout) <- c("sim", "id", "month", "therapy", "therapy_seq", 
                     "therapy_cycle", "death", "age", "ttd", "acr", "eular", 
                     "haq", "ttsi", "si", "yrlen")
  simout[, sim := sim + 1]
  simout[, id := id + 1]
  simout[, therapy_seq := therapy_seq + 1]
  simout[, therapy_cycle := therapy_cycle + 1]
  simout[, therapy := therapy + 1]
  return(simout)
}

#' Check parameters of sim_haq
#'
#' Error messages when incorrect inputs are passed to sim_haq.
#' 
#' @param input_data \code{input_data} as passed to \link{sim_haq}.
#' @param pars \code{pars} as passed to \link{sim_haq}.
#' @keywords internal
check_sim_haq <- function(input_data, pars){
  names.dist <- c("exponential", "exp", "weibull", "gompertz", "gamma", "llogis",
                  "lnorm", "gengamma")
  
  ## check input data
  if(is.null(input_data$haq0)) stop("'haq0' element of input_data list not given")
  if(is.null(input_data$age)) stop("'age' element of input_data list not given")
  if(is.null(input_data$male)) stop("'male' element of input_data list not given")
  if(is.null(input_data$x.mort)) stop("'x.mort' element of input_data list not given")
  if(is.null(input_data$x.dur)) stop("'x.dur' element of input_data list not given")
  n <- input_data$n
  if(length(input_data$age) != n | length(input_data$dis.dur) != n |
     nrow(input_data$x.mort) !=n | nrow(input_data$x.dur) != n) {
      stop(paste0("Number of patients not consistent accross elements of the input_data list.",
                  " Should equal ", n))
  }
  
  ## check parameter inputs
  n <- pars$n
  
  # rebound
  if(is.null(pars$rebound)) stop("'rebound' element of pars list not given")
  if(!is.vector(pars$rebound))  stop("'rebound' element of pars list must be a vector")
  if(length(pars$rebound) != n) {
    stop(paste0("Length of 'rebound' element of pars must be equal to the number of ",
         "sampled parameter sets which equals ", n))
  }
  
  # acr2eular
  if(is.null(pars$acr2eular)) stop("'acr2eular' element of pars list not given")
  if(!is.array(pars$acr2eular))  stop("'acr2eular' element of pars list must be an array")
  if(dim(pars$acr2eular)[1] != 4) stop("First dimension of 'acr2eular' element of pars must be 4")
  if(dim(pars$acr2eular)[2] != 3) stop("Second dimension of 'acr2eular' element of pars must be 3")
  if(dim(pars$acr2eular)[3] != n) {
    stop(paste0("Third dimension of 'acr2eular' element of pars must be equal to the number of",
         "sampled parameter sets which equals ", n))
  }

  # eular2haq
  if(is.null(pars$eular2haq)) stop("'eular2haq' element of pars list not given")
  if(!is.matrix(pars$eular2haq))  stop("'eular2haq' element of pars list must be a matrix")
  if(ncol(pars$eular2haq) != 3) stop("Number of columns of 'eular2haq' element of pars must be 3")
  if(nrow(pars$eular2haq) != n){
    stop(paste0("Number of rows in 'eular2haq' element of pars must be equal to the number ",
         "of sampled parameter sets which equals ", n))
  } 
  
  # haq.lprog.therapy
  if(is.null(pars$haq.lprog.therapy)) stop("'haq.lprog.therapy' element of pars list not given")
  if(!is.matrix(pars$haq.lprog.therapy))  stop("'haq.lprog.therapy' element of pars list must be a matrix")
  if(nrow(pars$haq.lprog.therapy) != n) {
    stop(paste0("Number of rows in 'haq.lprog.therapy' element of pars must be equal to the number ",
                "of sampled parameter sets which equals ", n))
  } 
  
  # haq.lprog.age
  if(is.null(pars$haq.lprog.age)) stop("'haq.lprog.age' element of pars list not given")
  if(!is.matrix(pars$haq.lprog.age))  stop("'haq.lprog.age' element of pars list must be a matrix")
  if(ncol(pars$haq.lprog.age) != 3) stop("Number of columns in 'haq.lprog.age' element of pars must be 3")
  if(nrow(pars$haq.lprog.age) != n) {
    stop(paste0("Number of rows in 'haq.lprog.age' element of pars must be equal to the number ",
                "of sampled parameter sets which equals ", n))
  }

  # logor.mort
  if(is.null(pars$logor.mort)) stop("'logor.mort' element of pars list not given")
  if(!is.matrix(pars$logor.mort))  stop("'logor.mort' element of pars list must be a matrix")
  if(nrow(pars$logor.mort) != n) {
    stop(paste0("Number of rows in 'logor.mort' element of pars must be equal to the number ",
                "of sampled parameter sets which equals ", n))
  }
  
  # mort.loghr.haqdif
  if(is.null(pars$mort.loghr.haqdif)) stop("'mort.loghr.haqdif' element of pars list not given")
  if(!is.matrix(pars$mort.loghr.haqdif))  stop("'mort.loghr.haqdif' element of pars list must be a matrix")
  if(ncol(pars$mort.loghr.haqdif) != 5) stop("Number of columns in 'mort.loghr.haqdif' element of pars must be 5")
  if(nrow(pars$mort.loghr.haqdif) != n) {
    stop(paste0("Number of rows in 'mort.loghr.haqdif' element of pars must be equal to the number ",
                "of sampled parameter sets which equals ", n))
  }
  
  # si.surv
  if(is.null(pars$si.surv)) stop("'si.surv' element of pars list not given")
  if(!is.list(pars$si.surv))  stop("'si.surv' element of pars list must be a list")
  if(all(!names(pars$si.surv) %in% names.dist)) {
    stop(paste0("Survival distribution in 'si.surv' element of pars list must contain at least ",
                "one of the following distributions: "), paste(names.dist, collapse = ", "))
  } 
  
  # ttd.eular.mod
  if(is.null(pars$ttd.eular.mod)) stop("'ttd.eular.mod' element of pars list not given")
  if(!is.list(pars$ttd.eular.mod))  stop("'ttd.eular.mod' element of pars list must be a list")
  if(all(!names(pars$ttd.eular.mod) %in% names.dist)) {
    dists.bad <- names(pars$ttd.eular.mod)[which(!names(pars$ttd.eular.mod) %in% names.dist)]
    stop(paste0("sim_haq does not support at least 1 of the survival distributions (",
                dists.bad,
                ") that are contained in the 'ttd.eular.mod'",
                " element of pars"))
  }  
  
  # ttd.eular.good
  if(is.null(pars$ttd.eular.good)) stop("'ttd.eular.good' element of pars list not given")
  if(!is.list(pars$ttd.eular.good))  stop("'ttd.eular.good' element of pars list must be a list")
  if(all(!names(pars$ttd.eular.good) %in% names.dist)) {
    dists.bad <- names(pars$ttd.eular.good)[which(!names(pars$ttd.eular.good) %in% names.dist)]
    stop(paste0("sim_haq does not support at least 1 of the survival distributions (",
                dists.bad,
                ") that are contained in the 'ttd.eular.good'",
                " element of pars"))
  }  
  
  # lt 
  if(is.null(pars$lt)) stop("'lt' element of pars list not given")
  if (!all(names(pars$lt) %in% c("female", "male"))) stop(paste0("'lt' element of pars must contain ",
                                                          "lifetables named 'male' and 'female'"))
  if(ncol(pars$lt$male) != 3 | ncol(pars$lt$female) != 3){
    stop("Number of columns in 'lt$female' element of pars and 'lt$male' must be 3")
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
#' model. See section 'pars' below.
#' @param check Should the function check parameters and input data passed to the model? Default is TRUE.
#' 
#' @section pars:
#' Parameters in \code{pars} are:
#' \describe{
#' \item{pain.mean}{Mean of pain score in the population.}
#' \item{haq.mean}{Mean of HAQ score in the population.}
#' \item{pain.var}{Variance of pain score in the population.}
#' \item{haq.var}{Variance of HAQ in the population.}
#' \item{painhaq.cor}{Correlation between pain and HAQ in the population.}
#' \item{mixture.utility}{Parameters sampled from the Hernandez Alva (2013) mixture model. See
#' the documentation in \link{sample_pars} for details.}
#' }
#' @return Matrix. First column is simulated pain score and second column is simulated utility. 
#' Each row corresponds to a unique patient and time-period (i.e. month) from \link{sim_haq}.
#'
#' @export
sim_utility_mixture <- function(simhaq, male, pars, check = TRUE){
  if (check) check_sim_utility_mixture(simhaq, male, pars)
  util <- sim_utility_mixtureC(simhaq$id - 1, simhaq$sim - 1, simhaq$haq, 
                   pars$pain.mean, pars$haq.mean, pars$pain.var, pars$haq.var, pars$painhaq.cor, 
                   simhaq$age, male,
                   pars$beta1, pars$beta2, pars$beta3, pars$beta4, 
                   pars$alpha1, pars$alpha2, pars$alpha3, pars$alpha4,
                   pars$alpha,
                   pars$epsilon1, pars$epsilon2, pars$epsilon3, pars$epsilon4,
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
    if(is.null(pars[[i]]))  stop(paste0("'", i, "'", " element of pars list not given"))
    if(!is.vector(pars$pain.mean)) stop(paste0("'", i, "'", " element of pars must be a vector"))
    if(length(pars$pain.mean) > 1) stop(paste0("'", i, "'", " element of pars must be length 1"))
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
#' @param input_data Variables other than those in simhaq needed to simulate model; these are
#' \code{Disease duration}, \code{Baseline HAQ}, \code{Male}, and \code{Previous DMARDs}. All variables
#' are generated from \link{input_data}.
#' @param pars Matrix of coefficients needed to simulate utility using the Wailoo (2006) model. 
#' See the documentation in \link{sample_pars} for details. Note that the matrix columns must be in 
#' the same order as generated by \link{sample_pars}. 
#' @param check Should the function check parameters and input data passed to the model? Default is TRUE.
#' 
#' @return Vector of utility scores. 
#'
#' @export
sim_utility_wailoo <- function(simhaq, input_data, pars, check = TRUE){
  if (check) check_sim_utility_wailoo(simhaq, input_data, pars)
  util <- sim_utility_wailooC(simhaq$sim - 1, simhaq$id - 1, simhaq$age,
                              input_data$dis.dur, input_data$haq0, input_data$male, 
                              input_data$prev.dmards, simhaq$haq,
                              pars[, "int"], pars[, "age"], pars[, "dis_dur"], 
                              pars[, "haq0"], pars[, "male"],
                              pars[, "prev_dmards"], pars[, "haq"])
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
check_sim_utility_wailoo <- function(simhaq, input_data, pars){
  # simhaq
  if(is.null(simhaq$id)) stop("'id' column of 'simhaq' is missing")
  if(is.null(simhaq$sim)) stop("'sim' column of 'simhaq' is missing")
  if(is.null(simhaq$age)) stop("'age' column of 'simhaq' is missing")
  if(is.null(simhaq$haq)) stop("'haq' column of 'simhaq' is missing")
  
  # input_data
  if(is.null(input_data$dis.dur)) stop("'dis.dur' element of input_data list not given")
  if(is.null(input_data$haq0)) stop("'haq0' element of input_data list not given")
  if(is.null(input_data$male)) stop("'male' element of input_data list not given")
  if(is.null(input_data$prev.dmards)) stop("'prev.dmards' element of input_data list not given")
  n <- input_data$n
  if(length(input_data$dis.dur) != n | length(input_data$haq0) != n |
     length(input_data$male) != n | length(input_data$prev.dmards) != n) {
    stop(paste0("Number of patients not consistent accross elements of the input_data list.",
                "Should equal ", n))
  }
  
  # pars
  n <- max(simhaq$sim)
  if(ncol(pars) != 7) stop("Number of columns in 'pars' must be 7")
  if(nrow(pars) != n) stop(paste0("Number of rows in 'pars' must be equal to the number of ",
                           "sampled parameter sets which equals ", n))
  if(any(!colnames(pars) %in% c("int", "age", "dis_dur", "haq0", "male", "prev_dmards", "haq"))){
    stop(paste0("Column names of 'pars' must be (in order) 'int', 'age', 'dis_dur'",
                "'haq0', 'male', 'prev_dmards', and 'haq'."))
  }
}

#' Simulate health care sector costs from simulated HAQ score
#'
#' Simulate health care sector costs from HAQ score simulated using \code{sim_haq}
#' 
#' @param simhaq Simulation output from \link{sim_haq}.
#' @param weight Weight of patient. Used to calculate treatment costs when dosing is based
#' on weight. 
#' @param pars List of sampled parameter values. The format is the same as the ouptut produced 
#' with \link{sample_pars}. The list must contain:
#' \describe{
#' \item{treat.cost}{Matrix of treatment cost data.}
#' \item{hosp.cost}{List of two matrices used to simulate hospital costs. }
#' \item{mgmt.cost}{Matrix of general management costs.}
#' }
#' @param cdmards.ind Index for cDMARDs. 
#' @param tcz.ind Index for tocilizumab.
#' @param tczmtx.ind Index for tocilizuma + methotrexate.
#' @param check Should the function check parameters and input data passed to the model? Default is TRUE.
#' 
#' @return Matrix. First column is simulated pain score and second column is simulated utility.
#'
#' @export
sim_hc_cost <- function(simhaq, weight, pars, 
                        cdmards.ind = which(therapy.pars$info[["sname"]] == "cdmards"),
                        tcz.ind = which(therapy.pars$info[["sname"]] == "tcz"),
                        tczmtx.ind = which(therapy.pars$info[["sname"]] == "tczmtx"),
                        check = TRUE){
  if (check) check_sim_hc_cost(simhaq, weight, pars)
  cycle_length <- 6
  treat.cost <- treat_costC(simhaq$therapy - 1, simhaq$therapy_cycle - 1, simhaq$id - 1, weight, cycle_length,
                            pars$treat.cost$ann_infusion_cost, pars$treat.cost$ann_rx_cost, 
                   pars$treat.cost$init_infusion_cost, pars$treat.cost$init_rx_cost, 
                   pars$treat.cost$weight_based, 
                   pars$treat.cost$ann_wgt_slope, pars$treat.cost$init_wgt_slope,
                   pars$treat.cost$ann_util, pars$treat.cost$init_util, 
                   pars$treat.cost$strength, pars$treat.cost$price, 
                   cdmards.ind - 1, tcz.ind - 1, tczmtx.ind - 1)
  treat.cost <- as.data.table(treat.cost)
  setnames(treat.cost, colnames(treat.cost), c("infusion_cost", "rx_cost"))
  treat.cost[, treat_cost := infusion_cost + rx_cost]
  hosp.cost <- hosp_costC(simhaq$haq, simhaq$yrlen, simhaq$sim - 1, 
                          pars$hosp.cost$cost.pday, pars$hosp.cost$hosp.days)
  mgmt.cost <- mgmt_costC(simhaq$yrlen, simhaq$sim - 1, rowSums(pars$mgmt.cost))
  si.cost <- si_costC(simhaq$si, simhaq$yrlen, simhaq$sim - 1, pars$si.cost)
  cost <- data.table(treat.cost, hosp_cost = hosp.cost, mgmt_cost = mgmt.cost,
                      si_cost = si.cost)
  cost[, hc_cost := treat_cost + hosp_cost + mgmt_cost + si_cost]
  return(cost)
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
  if(is.null(simhaq$therapy)) stop("'therapy' column of 'simhaq' is missing")
  if(is.null(simhaq$therapy_cycle)) stop("'therapy_cycle' column of 'simhaq' is missing")
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

#' Simulate productivity loss from simulated HAQ score
#'
#' Simulate productivity loss from HAQ score simulated using \code{sim_haq}
#' 
#' @param simhaq Simulation output from \link{sim_haq}.
#' @param pl_haq Decrease in wages (i.e. productivity loss) per one-unit increase in HAQ. Equivalent to output \code{prod.loss} from
#' \link{sample_pars}.
#' @return Vector of simulated productivity loss for each simulated patient and time-period.
#'
#' @export
sim_prod_loss <- function(simhaq, pl_haq, check = TRUE){
  if(check) check_prod_loss(simhaq, pl_haq)
  pl <- prod_lossC(simhaq$haq, simhaq$yrlen, simhaq$sim - 1, pl_haq)
  return(pl)
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


