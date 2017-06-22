rm(list = ls())
library("data.table")
mort.or <- fread("mort-or.csv")
mort.or <- mort.or[, .(var, or, or_se)] 
mort.or[, logor := log(or)]
mort.or[, logor_se := sqrt(or_se^2/or^2)] # delta method: var(exp(beta)) = exp(beta)var(beta)exp(beta)
mort.or[, logor_lower := logor - logor_se* qnorm(.975)]
mort.or[, logor_upper := logor + logor_se* qnorm(.975)]
save(mort.or, file = "../data/mort-or.rda", compress = "bzip2")