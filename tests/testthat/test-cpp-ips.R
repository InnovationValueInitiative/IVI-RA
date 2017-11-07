context("ips.cpp unit tests")
library("flexsurv")
library("data.table")
library("Rcpp")
seed <- runif(1, 0, 1000)

# Network-meta analyses -------------------------------------------------------
test_that("sim_acr_test", { 
  expect_equal(iviRA:::sim_acr_test(), 0)
})
test_that("sim_lm_test", { 
  expect_equal(iviRA:::sim_lm_test(), 0)
})

# Time to treatment discontinuation -------------------------------------------
test_that("ttd_da_test", { 
  x.da <- iviRA:::ttd_da_test(da_cat = 2)
  expect_equal(x.da[2], 1)
  x.da <- iviRA:::ttd_da_test(da_cat = 3)
  expect_equal(x.da[3], 1)
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
  expect_false(any(haq <0  | haq > 3))
})


# Test sim_tx_cost1C ----------------------------------------------------------
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

