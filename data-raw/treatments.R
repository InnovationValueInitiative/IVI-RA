rm(list = ls())
library(data.table)
treatments <- fread("treatments.csv")
save(treatments,
     file = "../data/treatments.rda", compress = "bzip2")

