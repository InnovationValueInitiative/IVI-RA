rm(list = ls())
library("data.table")
lifetable <- fread("lifetable.csv")
lifetable.male <- lifetable[male == 1, ]
lifetable.female <- lifetable[male == 0, ]
lifetable.male[, male := NULL]
lifetable.female[, male := NULL]
save(lifetable.male, file = "../data/lifetable-male.rda", compress = "bzip2")
save(lifetable.female, file = "../data/lifetable-female.rda", compress = "bzip2")