context("Workhorse IPS functions")
library("flexsurv")
library("data.table")
source("../../data-raw/func.R")

# Test C++ function durationC -------------------------------------------------
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
durationC_samp <- function(dist, type = 1){
  seed <- runif(1, 0, 100)
  fit <- fits[[dist]]
  est <- pars[[dist]]$est
  loc.est <- est[pars[[dist]]$loc.index]
  anc1.est <- est[pars[[dist]]$anc1.index]
  anc2.est <-  est[pars[[dist]]$anc2.index]

  set.seed(seed)
  samp1 <- IVI026::durationC(x, loc.est, anc1.est, loc.est, anc1.est, 
                             type = type, dist, cycle.length, 20,
                             anc2.est, anc2.est)
  set.seed(seed)
  samp2 <- IVI026::rsurvC(x %*% loc.est, anc1.est, dist, anc2.est)/cycle.length
  return(list(samp1, samp2))
}

test_that("durationC", {
  
  # exponential distribution
  samp <- durationC_samp("exp")
  expect_equal(samp[[1]], samp[[2]])
  samp <- durationC_samp("exp", type = 0)
  expect_equal(samp[[1]], 0)

  # weibull
  samp <- durationC_samp("weibull")
  expect_equal(samp[[1]], samp[[2]])
  
  # gompertz
  samp <- durationC_samp("gompertz")
  expect_equal(samp[[1]], samp[[2]])
  
  # log-normal
  samp <- durationC_samp("lnorm")
  expect_equal(samp[[1]], samp[[2]])
  
  # gamma
  samp <- durationC_samp("gamma")
  expect_equal(samp[[1]], samp[[2]])
  
  # log-logistic
  samp <- durationC_samp("llogis")
  expect_equal(samp[[1]], samp[[2]])
  
  # generalized gamma
  samp <- durationC_samp("gengamma")
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
  input.dat <- list(dis.dur = dis.dur, haq0 = haq0, male = male, prev.dmards = prev.dmards)
  
  x <- cbind(1, age, dis.dur, haq0, male, prev.dmards, haq)
  expect_equal(as.numeric(1/(1 + exp(-x %*% util.wailoo.pars$coef))),
              sim_utility_wailoo(simhaq, input.dat, t(as.matrix(util.wailoo.pars$coef))))
})


# small integration test ------------------------------------------------------
arms <- c(5, 14)
pat <- sample_pats(n = 1)
input.dat <- input_data(patdata = pat)
parsamp <- sample_pars(n = 100)
parsamp.table <- par_table(parsamp, pat)
sim.out <- sim_haq(arms, input_data = input.dat, pars = parsamp)
sim.out <- cbind(sim.out, sim_hc_cost(sim.out, pat[, "weight"], pars = parsamp))
sim.out[, prod_loss := sim_prod_loss(sim.out, pl_haq = parsamp$prod.loss)]
sim.out <- cbind(sim.out, sim_utility_mixture(sim.out, male = input.dat$male, 
                                      pars = c(pain, parsamp$mixture.utility)))
sim.out[, qalys := sim_qalys(sim.out, sim.out$utility, si_ul = parsamp$si.ul)]
util.wailoo <- sim_utility_wailoo(sim.out, input.dat, parsamp$wailoo.utility)
