rm(list = ls())
library("readxl")
library("data.table")
source("func.R")

# ACQUISITION AND ADMINISTRATION COSTS -----------------------------------------
treat.cost <- data.table(read_excel("cost.xlsx", sheet = "Cost"))
treat.lookup <- data.table(read_excel("cost.xlsx", sheet = "Lookup"))
treat.cost <- list(cost = treat.cost, lookup = treat.lookup)

# HOSPITALIZATION COSTS --------------------------------------------------------
haq <- c("0 to <0.5", "0.5 to <1", "1 to <1.5", "1.5 to <2", "2 to <2.5", ">2.5")
hosp.cost <- data.table(haq = haq,
                        days_mean = c(0.26, 0.13, 0.51, 0.72, 1.86, 4.16),
                        days_se = rep(0.5, 6),
                        cost_pday_mean = rep(1251, 6),
                        cost_pday_se = rep(191, 6))

# GENERAL MANAGEMENT COSTS -----------------------------------------------------
mgmt.cost <- fread("mgmt-cost.csv")
mgmt.cost[, se := se_normal(lower, upper, .975)]

# PRODUCTIVITY LOSS ------------------------------------------------------------
cpi <- data.table(read_excel("cpi.xlsx", range = "A12:O28"))
cpi.adj <- cpi[Year == 2017, HALF1]/cpi[Year == 2002, HALF1]
prod.loss.est <- 4372 * cpi.adj
prod.loss.se <- se_normal(2078 * cpi.adj, 6607 * cpi.adj)
prod.loss <- data.table(est = prod.loss.est, se = prod.loss.se)

# SAVE -------------------------------------------------------------------------
save(treat.cost, file = "../data/treat-cost.rda", compress = "bzip2")
save(hosp.cost, file = "../data/hosp-cost.rda", compress = "bzip2")
save(mgmt.cost, file = "../data/mgmt-cost.rda", compress = "bzip2")
save(prod.loss, file = "../data/prod-loss.rda", compress = "bzip2") 