rm(list = ls())
library("data.table")
source("func.R")
treatments <- fread("treatments.csv")

# SERIOUS INFECTION RATE -------------------------------------------------------
# Parameters from Cochrane review (Singh 2011)
ttsi <- data.table(sname = treatments$sname)
ttsi[, lograte := ifelse(sname %in% c("cdmards", "tt"), log(.026), 
                       ifelse(sname %in% c("placebo", "nbt"), log(0 + 10e-11),
                              log(.035)))]
ttsi[, lograte_se := se_normal(log(.027), log(.046), .975)]

# SAVE PARAMETERS --------------------------------------------------------------
save(ttsi,
     file = "../data/ttsi.rda", compress = "bzip2")
