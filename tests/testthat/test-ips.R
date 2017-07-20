context("Workhorse IPS functions")
library("flexsurv")
library("data.table")
source("../../data-raw/func.R")
seed <- runif(1, 0, 1000)

# select_model_structure ------------------------------------------------------
test_that("sim_duration_eular", {
  
  expect_error(
    select_model_structure(itreat_haq = "acr-haq", 
                           itreat_switch = "acr-eular-switch")
  )
  expect_error(
    select_model_structure(itreat_haq = "haq", 
                           itreat_switch = "acr-switch")
  )
  expect_error(
    select_model_structure(itreat_haq = "haq", 
                           itreat_switch = "acr-das28-switch")
  )
  expect_error(
    select_model_structure(itreat_haq = "haq", 
                           itreat_switch = "cr-sdai-switch")
  )
  expect_error(
    select_model_structure(itreat_haq = "haq", 
                           itreat_switch = "cr-sdai-switch")
  )
  expect_error(
    select_model_structure(itreat_haq = "haq", 
                           itreat_switch = "acr-eular-switch")
  )
  
})

# Test C++ function sim_duration ----------------------------------------------
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
test_sim_duration_eular <- function(dist, type = 1){
  fit <- fits[[dist]]
  est <- pars[[dist]]$est
  loc.est <- est[pars[[dist]]$loc.index]
  anc1.est <- est[pars[[dist]]$anc1.index]
  anc2.est <-  est[pars[[dist]]$anc2.index]

  set.seed(seed)
  samp1 <- iviRA::sim_duration_eular(x, loc.est, anc1.est, loc.est, anc1.est, 
                             type, dist, cycle.length, 20,
                             anc2.est, anc2.est)
  set.seed(seed)
  samp2 <- iviRA::rsurvC(x %*% loc.est, anc1.est, dist, anc2.est)/cycle.length
  return(list(samp1, samp2))
}

test_that("sim_duration_eular", {
  
  # exponential distribution
  samp <- test_sim_duration_eular("exp")
  expect_equal(samp[[1]], samp[[2]])
  samp <- test_sim_duration_eular("exp", type = 0)
  expect_equal(samp[[1]], 0)

  # weibull
  samp <- test_sim_duration_eular("weibull")
  expect_equal(samp[[1]], samp[[2]])
  
  # gompertz
  samp <- test_sim_duration_eular("gompertz")
  expect_equal(samp[[1]], samp[[2]])
  
  # log-normal
  samp <- test_sim_duration_eular("lnorm")
  expect_equal(samp[[1]], samp[[2]])
  
  # gamma
  samp <- test_sim_duration_eular("gamma")
  expect_equal(samp[[1]], samp[[2]])
  
  # log-logistic
  samp <- test_sim_duration_eular("llogis")
  expect_equal(samp[[1]], samp[[2]])
  
  # generalized gamma
  samp <- test_sim_duration_eular("gengamma")
  expect_equal(samp[[1]], samp[[2]])
})

# Test sim_utility_wailoo -----------------------------------------------------
test_that("sim_utility_wailoo", {
  age <- 55
  dis.dur <- 18.65
  haq0 <- 1
  male <- 1
  prev.dmards <- .0249
  haq <- 1.5
  simhaq <- data.table::data.table(sim = 1, id = 1, age = age, haq = haq)
  input.dat <- list(dis.dur = dis.dur, haq0 = haq0, male = male, prev.dmards = prev.dmards,
                    n = 1)
  x <- cbind(1, age, dis.dur, haq0, male, prev.dmards, haq)
  expect_equal(as.numeric(1/(1 + exp(-x %*% util.wailoo.pars$coef))),
              sim_utility_wailoo(simhaq, input.dat, t(as.matrix(util.wailoo.pars$coef))))
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
  beta.lc <- haq.lcgm.pars$coef[parameter == "beta3", coef]
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
  beta <- haq.lcgm.pars$coef[parameter %in% c("beta1", "beta2", "beta3", "beta4"), coef]
  beta <- matrix(beta, nrow = 4, ncol = 4, byrow = TRUE)
  age <- 55; female <- 1; das28 <- 6
  delta <- haq.lcgm.pars$coef[parameter %in% c("delta2", "delta3", "delta4"), coef] 
  delta <- matrix(delta, nrow = 3, ncol = length(delta)/3, byrow = TRUE)
  dhaq.C <- rep(NA, length(year))
  for (i in 1:length(year)){
    dhaq.C[i] <- iviRA:::sim_dhaq_lcgm1C(year[i], cycle.length, age, female, das28, delta, beta)
  }
  haq <- cumsum(c(1.5, dhaq.C))
  plot(c(2, year), haq)
})

# Test sim_itreat_haq ---------------------------------------------------------
test_that("sim_itreat_haq", {
  parsamp <- sample_pars(n = 3)
  line <- 0; therapy <- 0; nbt <- therapy + 5
  nma.acr1 <- nma.acr2 <- parsamp$acr$p1[1,, therapy + 1]
  nma.dhaq1 <- nma.dhaq2 <- parsamp$haq$dy1[1, therapy + 1]
  sim.acr2eular <- parsamp$acr2eular[,,1]
  sim.acr2haq <- parsamp$acr2haq
  sim.eular2haq <- parsamp$eular2haq[1, ]
  pars <- list("acr-haq", line, therapy, nbt, nma.acr1, nma.acr2,
            nma.dhaq1, nma.dhaq2, sim.acr2eular, sim.acr2haq, sim.eular2haq)
  
  ## Treatment -> ACR -> HAQ
  set.seed(seed)
  sim <- do.call(getFromNamespace("test_itreat_haq", "iviRA"), pars)
  set.seed(seed)
  acr <- hesim::rcat(matrix(nma.acr1, nrow = 1)) - 1
  eular <- hesim::rcat(acr2eular[acr + 1,, drop = FALSE]) - 1
  expect_equal(sim$acr, c(acr))
  expect_equal(sim$eular, c(eular))
  
  # check nbt returns correctly
  pars[[4]] <- therapy
  sim <- do.call(getFromNamespace("test_itreat_haq", "iviRA"), pars)
  expect_equal(sim$acr, 0)
  expect_equal(sim$eular, 0)
  expect_equal(sim$dhaq, 0)
  
  ## Treatment -> ACR -> EULAR -> HAQ
  pars[[1]] <- "acr-eular-haq"
  pars[[4]] <- nbt
  set.seed(seed)
  sim <- do.call(getFromNamespace("test_itreat_haq", "iviRA"), pars)
  set.seed(seed)
  acr <- c(hesim::rcat(matrix(nma.acr1, nrow = 1))) - 1
  eular <- c(hesim::rcat(acr2eular[acr + 1,, drop = FALSE])) - 1
  dhaq <- sim.eular2haq[eular + 1]
  expect_equal(sim$acr, acr)
  expect_equal(sim$eular, eular)
  expect_equal(sim$dhaq, as.vector(dhaq))
  
  # check nbt returns correctly
  pars[[4]] <- therapy
  sim <- do.call(getFromNamespace("test_itreat_haq", "iviRA"), pars)
  expect_equal(sim$acr, 0)
  expect_equal(sim$eular, 0)
  expect_equal(sim$dhaq, 0)
  
  ## Treatment -> HAQ
  pars[[1]] <- "haq"
  pars[[4]] <- nbt
  set.seed(seed)
  sim <- do.call(getFromNamespace("test_itreat_haq", "iviRA"), pars)
  set.seed(seed)
  expect_equal(sim$dhaq, as.vector(nma.dhaq1))
  expect_equal(sim$eular, eular)
  
  # check nbt returns correctly
  pars[[4]] <- therapy
  sim <- do.call(getFromNamespace("test_itreat_haq", "iviRA"), pars)
  expect_equal(sim$acr, 0)
  expect_equal(sim$eular, 0)
  expect_equal(sim$dhaq, 0)
})

# Test itreat_switchC ---------------------------------------------------------
# get_da_cat
test_that("get_da_cat", {
  expect_equal(iviRA:::get_das28_cat(2.1), 0)
  expect_equal(iviRA:::get_das28_cat(6.0), 3)
  expect_equal(iviRA:::get_sdai_cat(4.1), 1)
  expect_equal(iviRA:::get_sdai_cat(15.0), 2)
})

# itreat_switchC
test_that("itreat_switch", {
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
  sim <- do.call(getFromNamespace("test_itreat_switch", "iviRA"), pars)
  expect_equal(sim$t, 1)
  pars[["acr"]] <- 1
  sim <- do.call(getFromNamespace("test_itreat_switch", "iviRA"), pars)
  expect_equal(sim$tswitch, 0)
  
  # Treatment -> ACR -> DA -> Switch
  pars[[1]] <- "acr-das28-switch"
  set.seed(seed)
  sim <- do.call(getFromNamespace("test_itreat_switch", "iviRA"), pars)
  set.seed(seed)
  twitch.R <- rbinom(1, 1, p[iviRA:::get_das28_cat(das28 + sim.acr2das28[pars[["acr"]] + 1])])
  expect_equal(sim$tswitch, twitch.R)
  
  # Treatment -> DA -> Switch
  pars[[1]] <- "das28-switch"
  set.seed(seed)
  sim <- do.call(getFromNamespace("test_itreat_switch", "iviRA"), pars)
  set.seed(seed)
  twitch.R <- rbinom(1, 1, p[iviRA:::get_das28_cat(das28 + nma.das28.1[therapy + 1])])
  expect_equal(sim$tswitch, twitch.R)
  

})

# small integration test ------------------------------------------------------
pop <- sample_pats(n = 10)
arm.names <- c("adamtx", "cdmards")
parsamp <- sample_pars(n = 100)
mod.struct <- select_model_structure(itreat_haq = "acr-haq",
                                     itreat_switch = "acr-switch",
                                     cdmards_haq_model = "lcgm",
                                     utility_model = "wailoo")
input.dat <- get_input_data(patdata = pop, model_structure = mod.struct)
sim.out <- sim_iviRA(arms = arm.names, input_data = input.dat, pars = parsamp,
                     model_structure = mod.struct)



