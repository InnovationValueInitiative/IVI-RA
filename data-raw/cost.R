rm(list = ls())
library("readxl")

# IMPORT COSTS -----------------------------------------------------------------
ra.cost <- data.frame(read_excel("data-raw/cost.xlsx", sheet = "RData"))
save(ra.cost, 
     file = "data/cost.rda", compress = "bzip2")