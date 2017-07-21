rm(list = ls())
library("readxl")
library("data.table")
source("func.R")

# ACQUISITION AND ADMINISTRATION COSTS -----------------------------------------
treat.cost <- data.table(read_excel("cost.xlsx", sheet = "RData"))

# GENERAL MANAGEMENT COSTS -----------------------------------------------------
mgmt.cost <- fread("mgmt-cost.csv")
mgmt.cost[, se := se_normal(lower, upper, .975)]

# SAVE -------------------------------------------------------------------------
save(treat.cost, file = "../data/treat-cost.rda", compress = "bzip2")
save(mgmt.cost, file = "../data/mgmt-cost.rda", compress = "bzip2")