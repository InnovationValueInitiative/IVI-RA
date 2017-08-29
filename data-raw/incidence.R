rm(list = ls())
library("data.table")

# load data
incidence.grps <- fread("incidence.csv")
incidence.grps[, person_years := events /(incidence_per_100000/100000)]
incidence.grps[, nrows := age_top - age_bot + 1]
incidence.grps <- incidence.grps[rep(seq_len(nrow(incidence.grps)), times = nrows)]
incidence.grps[, age := rep(seq(18, 85), times = 2)]
incidence.grps[, incidence_rate := incidence_per_100000/100000]
incidence.grps <- incidence.grps[, .(age, gender, events, person_years, incidence_rate)]

# female
incidence.female <- incidence.grps[gender == "female"]
incidence.female[, gender := NULL]
incidence.female <- rbind(data.table(age = seq(0, 17), events = 0, person_years = 0,
                                     incidence_rate = 0), incidence.female)

# male
incidence.male <- incidence.grps[gender == "male"]
incidence.male[, gender := NULL]
incidence.male <- rbind(data.table(age = seq(0, 17), events = 0, person_years = 0,
                                     incidence_rate = 0), incidence.male)

# save
save(incidence.female, file = "../data/incidence-female.rda", compress = "bzip2")
save(incidence.male, file = "../data/incidence-male.rda", compress = "bzip2")


