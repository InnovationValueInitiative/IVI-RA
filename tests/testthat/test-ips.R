context("Workhorse IPS functions")
library("flexsurv")
library("data.table")
library("Rcpp")
source("../../data-raw/func.R")
seed <- runif(1, 0, 1000)

# select_model_structure ------------------------------------------------------
test_that("select_model_structures", {
  
  # default options are constant when selected options are greater than length 1
  expect_error(
    select_model_structures(tx_ihaq = c("acr-eular-haq", "acr-haq"),
                            tx_iswitch = c("acr-eular-switch", "acr-das28-switch")), NA
  )
  
  expect_warning(
    select_model_structures(tx_ihaq = c("acr-eular-haq", "acr-haq"),
                            tx_iswitch = c("acr-eular-switch", "acr-das28-switch")), NA
  )
  
  # return error when options are of different length
  expect_error(
    select_model_structures(tx_ihaq = c("acr-eular-haq", "acr-haq", "haq"),
                            tx_iswitch = c("acr-eular-switch", "acr-das28-switch"))
  )
  
  # combinations of model structures must be correct
  expect_error(
    select_model_structures(tx_ihaq = c("acr-eular-haq", "acr-haq"),
                           tx_iswitch = c("acr-eular-switch", "acr-eular-switch"))
  )
  expect_error(
    select_model_structures(tx_ihaq = "haq", 
                           tx_iswitch = "acr-switch")
  )
  expect_error(
    select_model_structures(tx_ihaq = "haq", 
                           tx_iswitch = "acr-das28-switch")
  )
  expect_error(
    select_model_structures(tx_ihaq = "haq", 
                           tx_iswitch = "acr-sdai-switch")
  )
  expect_error(
    select_model_structures(tx_ihaq = "haq", 
                           tx_iswitch = "acr-sdai-switch")
  )
  expect_error(
    select_model_structures(tx_ihaq = "haq", 
                           tx_iswitch = "acr-eular-switch")
  )
  
})

# Test treatment duration -----------------------------------------------------
# Data for testing
dists <- c("exp", "weibull", "gompertz", "lnorm", "llogis", "gamma", "gengamma")
fits <- pars <- vector(mode = "list", length(dists))
names(pars) <- names(fits) <- dists
for (i in 1:length(dists)){
  fits[[i]] <- suppressWarnings(flexsurv::flexsurvreg(formula = Surv(futime, fustat) ~ 1 + age, 
                                                     data = ovarian, dist = dists[i]))
  pars[[i]] <- flexsurvreg_pars(fits[[i]])
}
x <- c(1, ovarian[1, "age"])
cycle.length <- 6

# test function
test_sim_ttd_eular <- function(dist, type = 1){
  fit <- fits[[dist]]
  est <- pars[[dist]]$est
  loc.est <- est[pars[[dist]]$loc.index]
  anc1.est <- est[pars[[dist]]$anc1.index]
  anc2.est <-  est[pars[[dist]]$anc2.index]

  set.seed(seed)
  samp1 <- iviRA::sim_ttd_eular(x, loc.est, anc1.est, loc.est, anc1.est, 
                             type, dist, cycle.length, 20,
                             anc2.est, anc2.est)
  set.seed(seed)
  samp2 <- iviRA::rsurvC(x %*% loc.est, anc1.est, dist, anc2.est)/cycle.length
  return(list(samp1, samp2))
}

test_that("sim_ttd_eular", {
  
  # exponential distribution
  samp <- test_sim_ttd_eular("exp")
  expect_equal(samp[[1]], samp[[2]])
  samp <- test_sim_ttd_eular("exp", type = 0)
  expect_equal(samp[[1]], 0)

  # weibull
  samp <- test_sim_ttd_eular("weibull")
  expect_equal(samp[[1]], samp[[2]])
  
  # gompertz
  samp <- test_sim_ttd_eular("gompertz")
  expect_equal(samp[[1]], samp[[2]])
  
  # log-normal
  samp <- test_sim_ttd_eular("lnorm")
  expect_equal(samp[[1]], samp[[2]])
  
  # gamma
  samp <- test_sim_ttd_eular("gamma")
  expect_equal(samp[[1]], samp[[2]])
  
  # log-logistic
  samp <- test_sim_ttd_eular("llogis")
  expect_equal(samp[[1]], samp[[2]])
  
  # generalized gamma
  samp <- test_sim_ttd_eular("gengamma")
  expect_equal(samp[[1]], samp[[2]])
})

# Test sim_dhaq_lcgm1C --------------------------------------------------------
test_that("sim_mlogit_classC", {
  delta <- matrix(c(seq(1, 6)), nrow = 2, ncol = 3)
  w <- c(1, 1, 2)
  set.seed(seed)
  latclass <- iviRA:::sim_mlogit_classC(w, delta)
  expect_is(latclass, "integer")
})

test_that("sim_haq_class_lcgm1C", {
  beta.lc <- iviRA::haq.lcgm$coef[parameter == "beta3", est]
  year <- seq(2.5, 10, .5)
  year.lag <- year - .5
  xt <- 1 - (1/(1 + year))
  xt.lag <-  1 - (1/(1 + year.lag))

  # compare R with C++
  dhaq.C <- rep(NA, length(xt))
  for (i in 1:length(xt)){
    dhaq.C[i] <- iviRA:::sim_dhaq_class_lcgm1C(year[i], 6, beta.lc)
  }
  dhaq.R <- beta.lc[2] * (xt - xt.lag) + beta.lc[3] * (xt^2 - xt.lag^2) +
    beta.lc[4] * (xt^3 - xt.lag^3)
  expect_equal(dhaq.R, dhaq.C)
  
  # HAQ over time
  haq <- cumsum(c(1.5, dhaq.C))
  plot(year, haq[-1])
})

test_that("sim_haq_lcgm1C", {
  year <- seq(2.5, 10); cycle.length <- 6
  beta <- iviRA::haq.lcgm$coef[parameter %in% c("beta1", "beta2", "beta3", "beta4"), est]
  beta <- matrix(beta, nrow = 4, ncol = 4, byrow = TRUE)
  age <- 55; female <- 1; das28 <- 6
  delta <- haq.lcgm$coef[parameter %in% c("delta2", "delta3", "delta4"), est] 
  delta <- matrix(delta, nrow = 3, ncol = length(delta)/3, byrow = TRUE)
  dhaq.C <- rep(NA, length(year))
  for (i in 1:length(year)){
    dhaq.C[i] <- iviRA:::sim_dhaq_lcgm1C(year[i], cycle.length, age, female, das28, delta, beta)
  }
  haq <- cumsum(c(1.5, dhaq.C))
  plot(c(2, year), haq)
})

# Test sim_tx_ihaq ---------------------------------------------------------
test_that("sim_tx_ihaq", {
  parsamp <- sample_pars(n = 3)
  line <- 0; therapy <- 0; nbt <- therapy + 5
  nma.acr1 <- nma.acr2 <- parsamp$acr$p1[1,, therapy + 1]
  nma.dhaq1 <- nma.dhaq2 <- parsamp$haq$dy1[1, therapy + 1]
  sim.acr2eular <- parsamp$acr2eular[,,1]
  sim.acr2haq <- parsamp$acr2haq
  sim.eular2haq <- parsamp$eular2haq[1, ]
  pars <- list("acr-haq", line, therapy, nbt, nma.acr1, nma.acr2,
            nma.dhaq1, nma.dhaq2, sim.acr2eular, sim.acr2haq, sim.eular2haq)
  
  # Treatment -> ACR -> HAQ
  set.seed(seed)
  sim <- do.call(getFromNamespace("test_tx_ihaq", "iviRA"), pars)
  set.seed(seed)
  acr <- hesim::rcat(matrix(nma.acr1, nrow = 1)) - 1
  expect_equal(sim$acr, c(acr))
  expect_true(is.na(sim$eular))
  
  ## check nbt returns correctly
  pars[[4]] <- therapy
  sim <- do.call(getFromNamespace("test_tx_ihaq", "iviRA"), pars)
  expect_equal(sim$acr, 0)
  expect_true(is.na(sim$eular))
  expect_equal(sim$dhaq, 0)
  
  # Treatment -> ACR -> EULAR -> HAQ
  pars[[1]] <- "acr-eular-haq"
  pars[[4]] <- nbt
  set.seed(seed)
  sim <- do.call(getFromNamespace("test_tx_ihaq", "iviRA"), pars)
  set.seed(seed)
  acr <- c(hesim::rcat(matrix(nma.acr1, nrow = 1))) - 1
  eular <- c(hesim::rcat(acr2eular[acr + 1,, drop = FALSE])) - 1
  dhaq <- sim.eular2haq[eular + 1]
  expect_equal(sim$acr, acr)
  expect_equal(sim$eular, eular)
  expect_equal(sim$dhaq, as.vector(dhaq))
  
  ## check nbt returns correctly
  pars[[4]] <- therapy
  sim <- do.call(getFromNamespace("test_tx_ihaq", "iviRA"), pars)
  expect_equal(sim$acr, 0)
  expect_equal(sim$eular, 0)
  expect_equal(sim$dhaq, 0)
  
  # Treatment -> HAQ
  pars[[1]] <- "haq"
  pars[[4]] <- nbt
  set.seed(seed)
  sim <- do.call(getFromNamespace("test_tx_ihaq", "iviRA"), pars)
  set.seed(seed)
  expect_true(is.na(sim$acr))
  expect_true(is.na(sim$eular))
  expect_equal(sim$dhaq, as.vector(nma.dhaq1))
  
  ## check nbt returns correctly
  pars[[4]] <- therapy
  sim <- do.call(getFromNamespace("test_tx_ihaq", "iviRA"), pars)
  expect_true(is.na(sim$acr))
  expect_true(is.na(sim$eular))
  expect_equal(sim$dhaq, 0)
})

# Test tx_iswitchC ---------------------------------------------------------
test_that("get_da_cat", {
  expect_equal(iviRA:::get_das28_cat(2.1), 0)
  expect_equal(iviRA:::get_das28_cat(6.0), 3)
  expect_equal(iviRA:::get_sdai_cat(4.1), 1)
  expect_equal(iviRA:::get_sdai_cat(15.0), 2)
})

test_that("tx_iswitch", {
  parsamp <- sample_pars(n = 3)
  line <- 0; therapy <- 1; nbt <- therapy + 2; 
  acr <- 0; eular <- 2; das28 <- 6; cdai <- 41; sdai <- 43
  sim.acr2das28 <- parsamp$acr2das28[1, ]
  sim.acr2sdai <- parsamp$acr2sdai[1, ]
  sim.acr2cdai <- parsamp$acr2cdai[1, ]
  nma.das28.1 <- parsamp$das28$dy1[1, ]
  nma.das28.2 <- parsamp$das28$dy2[1, ]
  p <- c(0, .15, .22, .29)
  pars <- list("acr-switch", line, therapy, nbt, acr = acr, eular, das28, sdai, cdai,
            sim.acr2das28, sim.acr2sdai, sim.acr2cdai, nma.das28.1, nma.das28.2, p)
  
  # Treatment -> ACR -> Switch
  sim <- do.call(getFromNamespace("test_tx_iswitch", "iviRA"), pars)
  expect_equal(sim$t, 1)
  pars[["acr"]] <- 1
  sim <- do.call(getFromNamespace("test_tx_iswitch", "iviRA"), pars)
  expect_equal(sim$tswitch, 0)
  
  # Treatment -> ACR -> DA -> Switch
  pars[[1]] <- "acr-das28-switch"
  set.seed(seed)
  sim <- do.call(getFromNamespace("test_tx_iswitch", "iviRA"), pars)
  set.seed(seed)
  twitch.R <- rbinom(1, 1, p[iviRA:::get_das28_cat(das28 + sim.acr2das28[pars[["acr"]] + 1])])
  expect_equal(sim$tswitch, twitch.R)
  
  # Treatment -> DA -> Switch
  pars[[1]] <- "das28-switch"
  set.seed(seed)
  sim <- do.call(getFromNamespace("test_tx_iswitch", "iviRA"), pars)
  set.seed(seed)
  twitch.R <- rbinom(1, 1, p[iviRA:::get_das28_cat(das28 + nma.das28.1[therapy + 1])])
  expect_equal(sim$tswitch, twitch.R)
})

# Test simulated costs --------------------------------------------------------
gen_agents <- function(name){
  tx.lookup <- iviRA::tx.cost$lookup[sname == name]
  agents <- matrix(match(unlist(tx.lookup[, -1, with = FALSE]), 
                            iviRA::tx.cost$cost$sname) - 1,
                      nrow = nrow(tx.lookup))
  return(agents)
}
agents <- gen_agents("etnmtx")
tc <- iviRA::tx.cost$cost
discount <- .25
pars.tc <- list(t = 0, agents_ind = agents, tx_name = tc$sname, 
             init_dose_val = tc$init_dose_val, ann_dose_val = tc$ann_dose_val,
             strength_val = tc$strength_val, 
             init_num_doses = tc$init_num_doses, ann_num_doses = tc$ann_num_doses,
             price = tc$price_per_unit, 
             infusion_cost = tc$infusion_cost, loading_dose = tc$loading_dose,
             weight_based = tc$weight_based, weight = 112, cycle_length = 6, 
             discount = rep(discount, nrow(tc)))

test_that("sim_tx_cost1C", {
  t <- 0
  weight <- 112
  
  # etnmtx
  ## first 6 months
  tcC <- do.call(getFromNamespace("sim_tx_cost1C", "iviRA"), pars.tc)
  tcR <- sum((tc[sname %in% c("etn", "cdmards"), .(cost = ceiling(init_dose_val/strength_val) * 
                                   init_num_doses * price_per_unit * (1 - discount))]))
  expect_equal(tcC, tcR)
  
  ## maintenance phase
  pars.tc$t <- 4
  tcC <- do.call(getFromNamespace("sim_tx_cost1C", "iviRA"), pars.tc)
  tcR <- sum((tc[sname %in% c("etn", "cdmards"), 
                 .(cost = ceiling(ann_dose_val/strength_val) * 
                   ann_num_doses * price_per_unit * (1 - discount))]))
  expect_equal(tcC, tcR/2)
  
  # tczmtx
  agents <- gen_agents("tczmtx")
  pars.tc$agents_ind <- agents
  
  ## over 100 kg
  tcC <- do.call(getFromNamespace("sim_tx_cost1C", "iviRA"), pars.tc)
  tcR <- as.numeric(tc[sname =="tcz", 
                       .(cost = ceiling(ann_dose_val/strength_val) * 
                                ann_num_doses * 2 * price_per_unit * (1 - discount))] +
    tc[sname =="cdmards", .(cost = ceiling(ann_dose_val/strength_val) * 
                          ann_num_doses  * price_per_unit * (1 - discount))])
  expect_equal(tcC, tcR/2)
  
  ## less than 100 kg
  pars.tc$weight <- 85
  tcC <- do.call(getFromNamespace("sim_tx_cost1C", "iviRA"), pars.tc)
  tcR <- as.numeric(tc[sname =="tcz", 
                       .(cost = ceiling(ann_dose_val/strength_val) * 
                                ann_num_doses * price_per_unit * (1 - discount))] +
                      tc[sname =="cdmards", 
                         .(cost = ceiling(ann_dose_val/strength_val) * 
                            ann_num_doses  * price_per_unit * (1 - discount))])
  expect_equal(tcC, tcR/2)
  
  # ifxmtx
  agents <- gen_agents("ifxmtx")
  pars.tc$agents_ind <- agents
  
  ## maintenance phase
  tcC <- do.call(getFromNamespace("sim_tx_cost1C", "iviRA"), pars.tc)
  tcR <- as.numeric(tc[sname =="ifx", 
                       .(cost = ceiling(pars.tc$weight * ann_dose_val/strength_val) * 
                                          ann_num_doses * price_per_unit * (1 - discount) +
                                          ann_num_doses * infusion_cost)] +
                      tc[sname =="cdmards", 
                         .(cost = ceiling(ann_dose_val/strength_val) * 
                           ann_num_doses  * price_per_unit * (1 - discount))])
  expect_equal(tcC, tcR/2)
  
  # abtscmtx
  agents <- gen_agents("abtscmtx")
  pars.tc$agents_ind <- agents
  
  ## maintenance phase
  tcC <- do.call(getFromNamespace("sim_tx_cost1C", "iviRA"), pars.tc)
  tcR <- sum((tc[sname %in% c("abtsc", "cdmards"), 
                 .(cost = ceiling(ann_dose_val/strength_val) * 
                  ann_num_doses * price_per_unit * (1 - discount))]))
  expect_equal(tcC, tcR/2)
})

# hospital, management, serious infection, and productivity costs 
haq <- 1.7
yrlen <- .5
parsamp <- sample_pars(n = 10)

test_that("sim_hosp_cost1C", {
  expect_equal(iviRA:::sim_hosp_cost1C(haq, yrlen, parsamp$hosp.cost$hosp.days[1, ], 
                                 parsamp$hosp.cost$cost.pday[1, ]),
  as.vector(parsamp$hosp.cost$hosp.days[1, 4] * parsamp$hosp.cost$cost.pday[1, 4] * 
              yrlen))
})

test_that("sim_mgmt_cost1C", {
  expect_equal(iviRA:::sim_mgmt_cost1C(yrlen, sum(parsamp$mgmt.cost[1, ])),
               sum(parsamp$mgmt.cost[1, ]) * yrlen)
})

test_that("sim_si_cost1C", {
  expect_equal(iviRA:::sim_si_cost1C(0, yrlen, parsamp$si.cost[1]),
               0)
  expect_equal(iviRA:::sim_si_cost1C(1, yrlen, parsamp$si.cost[1]),
               parsamp$si.cost[1])
})

test_that("sim_prod_loss1C", {
  expect_equal(iviRA:::sim_prod_loss1C(haq, yrlen, parsamp$prod.loss[1]),
              parsamp$prod.loss[1] * haq *  yrlen)
})

# Test utility simulation(s) --------------------------------------------------
# Wailoo (2006)
test_that("sim_utility_wailoo1C", {
  age <- 55; dis.dur <- 18.65; haq0 <- 1; male <- 1
  prev.dmards <- .0249; haq <- 1.5
  simhaq <- data.table::data.table(sim = 1, id = 1, age = age, haq = haq)
  x <- cbind(1, age, dis.dur, haq0, male, prev.dmards, haq)
  
  # c++ double
  util.C1 <- iviRA:::sim_utility_wailoo1C(age = age, disease_duration = dis.dur,
                                          haq0 = haq0, male = male, 
                                          prev_dmards = prev.dmards, haq = haq,
                                          b = iviRA::utility.wailoo$est)
  
  # c++ vector
  beta <- t(as.matrix(iviRA::utility.wailoo$est))
  colnames(beta) <- c("int", "age", "dis_dur", "haq0", "male", "prev_dmards", "haq")
  util.C <- sim_utility_wailoo(simhaq = simhaq, haq0 = haq0, male = male, 
                               prev_dmards = prev.dmards,
                               coefs = beta)
  
  # R
  util.R <- as.numeric(1/(1 + exp(-x %*% iviRA::utility.wailoo$est)))
  
  # Compare
  expect_equal(util.C1, util.R)
  expect_equal(util.C, util.R)
})

# Hernandez-Alava (2013)
# test_that("sim_utility_mixture1C", {
#   haq <- 1.5; age <- 55; male <- 1
#   pars <- sample_pars(n = 10)
#   pars <- pars$utility.mixture
#   iviRA:::sim_utility_mixture1C(haq = haq, pain.mean = pars$pain$pain.mean,
#                         haq_mean = pars$pain$haq.mean, pain_var = pars$pain$pain.var,
#                         haq_var = pars$pain$haq.var, painhaq_cor = pars$pain$painhaq.cor,
#                         age = age, male = male, 
#                         beta1 = pars$beta1[1, ], beta2 = pars$beta2[1, ], 
#                         beta3 = pars$beta3[1, ], beta4 = pars$beta4[1, ], 
#                         alpha1 = pars$alpha1[1], alpha2 = pars$alpha2[1],
#                         alpha3 = pars)
# })

# Test sim_qalys --------------------------------------------------------------
n <- 10
parsamp <- sample_pars(n = 5)
x.attr <- iviRA::tx.attr$data
simhaq <- data.table(yrlen = rep(.5, n), sim = rep(seq(1, 5), each = n/5),
                     tx = which(iviRA::treatments$sname == "tof"),
                     si = rbinom(n, 1, .1))
utility <- runif(n, 0, 1)
tx.attr.ug <- parsamp$tx.attr.ug 
tx.attr.ug[, 1] <- 0.1

test_that("sim_qalys", {
  sim.qalys <- sim_qalys(simhaq = simhaq, utility = utility,
                         si_ul = parsamp$si.ul, x_attr = x.attr, 
                         tx_attr_ug = tx.attr.ug)
  sim.qalys1 <- simhaq$yrlen[1] * (utility[1] - simhaq$si[1] * parsamp$si.ul[1]/12 +
                      x.attr[simhaq$tx[1],, drop = FALSE] %*% t(tx.attr.ug[1,, drop = FALSE]))
  expect_equal(sim.qalys[1], sim.qalys1[1])
  
  utility[1] <- 2
  sim.qalys <- sim_qalys(simhaq = simhaq, utility = utility,
                         si_ul = parsamp$si.ul, x_attr = x.attr, 
                         tx_attr_ug = tx.attr.ug)
  expect_equal(sim.qalys[1], 0.5)
})

# Summary output for individual means -----------------------------------------
# mod.SimMeans <- Rcpp::Module('mod_SimMeans', PACKAGE = "iviRA")
# SimMeans <- mod.SimMeans$SimMeans
# sm <- new(SimMeans, 2, 5, 10, .03, .03)
# sm$get_id()
# sm$get_varsums()
# sm$increment_id(0, 2, 1)
# sm$increment_varsums(2.0)
# sm$get_id()
# sm$get_varsums()
# sm$calc_means()

# Summary output for time means -----------------------------------------------
# mod.TimeMeans <- Rcpp::Module('mod_TimeMeans', PACKAGE = "iviRA")
# TimeMeans <- mod.TimeMeans$TimeMeans
# tm <- new(TimeMeans, 2, 3, 4, 5, 6.0)
# tm$get_id()
# tm$get_alive()
# tm$increment_id(1, 1, 2)
# tm$increment_alive()
# tm$get_alive()
# tm$increment_varsums(.3, .5)
# tm$get_id()
# tm$get_index()
# tm$get_varsums()
# tm$calc_means()

# Summary output during model cycle 0 -----------------------------------------
# mod.Out0 <- Rcpp::Module('mod_Out0', PACKAGE = "iviRA")
# Out0 <- mod.Out0$Out0
# out <- new(Out0, 2, 5, 10, 5)
# out$push_back(0, 0, 0, 0, 0, 1, 2, 2.4, 2.9)
# out$get_acr()
# out$get_ttsi()

# small integration test ------------------------------------------------------
pop <- sample_pats(n = 100)
arm.names <- c("adamtx", "cdmards")
mod.structs <- select_model_structures(tx_ihaq = c("acr-haq", "acr-eular-haq"),
                                     tx_iswitch = c("acr-switch", "acr-eular-switch"),
                                     cdmards_haq_model = c("lcgm", "linear"),
                                     ttd_dist = c("gengamma", "lnorm"),
                                     utility_model = c("mixture", "wailoo"))
input.dat <- get_input_data(patdata = pop, model_structure = mod.structs)
parsamp <- sample_pars(n = 100)
sim.out <- sim_iviRA(arms = arm.names, input_data = input.dat, pars = parsamp,
                     model_structures = mod.structs, output = "summary")


