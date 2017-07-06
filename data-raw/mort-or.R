rm(list = ls())
library("data.table")
source("func.R")
mort.or <- fread("mort-or.csv")
mort.or[, logor := log(or)]
mort.or[, logor_lower := log(or_lower)]
mort.or[, logor_upper := log(or_upper)]
mort.or[, logor_se := sqrt(or_se^2/or^2)] 
logor.se.dm <- sqrt(mort.or$or_se^2/mort.or$or^2) # delta method: var(exp(beta)) = exp(beta)var(beta)exp(beta)
mort.or[, logor_se := se_normal(logor_lower, logor_upper)]
save(mort.or, file = "../data/mort-or.rda", compress = "bzip2")