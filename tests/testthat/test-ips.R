context("Workhorse IPS functions")
library("flexsurv")
library("data.table")
library("Rcpp")
source("../../data-raw/func.R")

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

# Test sim_qalys --------------------------------------------------------------
n <- 10
pop <- sample_pop(n = 10)
input.dat <- get_input_data(pop)
parsamp <- sample_pars(n = 5, input_data = input.dat)
x.attr <- iviRA::utility.tx.attr$x
simhaq <- data.table(yrlen = rep(.5, n), sim = rep(seq(1, 5), each = n/5),
                     tx = which(iviRA::treatments$sname == "tof"),
                     si = rbinom(n, 1, .1))
utility <- runif(n, 0, 1)
tx.attr.coef <- parsamp$utility.tx.attr
tx.attr.coef[, 1] <- 0.1

test_that("sim_qalys", {
  sim.qalys <- sim_qalys(simhaq = simhaq, utility = utility,
                         si_ul = parsamp$si.ul, x_attr = x.attr, 
                         tx_attr_coef = tx.attr.coef)
  sim.qalys1 <- simhaq$yrlen[1] * (utility[1] - simhaq$si[1] * parsamp$si.ul[1]/12 +
                      x.attr[simhaq$tx[1],, drop = FALSE] %*% t(tx.attr.ug[1,, drop = FALSE]))
  expect_equal(sim.qalys[1], sim.qalys1[1])
  
  utility[1] <- 2
  sim.qalys <- sim_qalys(simhaq = simhaq, utility = utility,
                         si_ul = parsamp$si.ul, x_attr = x.attr, 
                         tx_attr_ug = tx.attr.ug)
  expect_equal(sim.qalys[1], 0.5)
})

# sim_iviRA -------------------------------------------------------------------
pop <- sample_pop(n = 10, type = "homog")
arm.names <- c("adamtx", "cdmards")
mod.structs <- select_model_structures(tx_ihaq = c("acr-haq", "acr-eular-haq"),
                                     tx_iswitch = c("acr-switch", "acr-eular-switch"),
                                     cdmards_haq_model = c("lcgm", "linear"),
                                     ttd_cause = c("all", "si"),
                                     ttd_dist = c("gengamma", "lnorm"),
                                     utility_model = c("mixture", "wailoo"))
input.dat <- get_input_data(pop = pop)
parsamp <- sample_pars(n = 10, input_dat = input.dat)
sim.out <- sim_iviRA(arms = arm.names, input_data = input.dat, pars = parsamp,
                     model_structures = mod.structs, output = "data")


