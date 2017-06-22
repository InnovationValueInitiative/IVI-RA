rm(list = ls())
library("data.table")
source("func.R")
mgmt.cost <- fread("mgmt-cost.csv")
mgmt.cost[, se := se_normal(lower, upper, .975)]
save(mgmt.cost, 
     file = "../data/mgmt-cost.rda", compress = "bzip2")