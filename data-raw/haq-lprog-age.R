rm(list = ls())
library("data.table")
source("func.R")

# results from table 1 in Michaud (2011)
overall.mean <- .014
overall.lower <- .012
overall.upper <- .015

agel40.mean <- -.006
agel40.lower <- -.014
agel40.upper <- .002

age40to65.mean <- .006
age40to65.lower <- .004
age40to65.upper <- .007

age65p.mean <- .031
age65p.lower <- .011
age65p.upper <- .019

# calculate standard errors 
overall.se <- se_normal(overall.lower, overall.upper)
agel40.se <- se_normal(agel40.lower, agel40.upper)
age40to65.se <- se_normal(age40to65.lower, age40to65.upper)
age65p.se <- se_normal(age65p.lower, age65p.upper)

# differences
agel40.diff <- agel40.mean - overall.mean
agel40.diff.se <- sqrt(agel40.se^2 + overall.se^2)
age40to65.diff <- age40to65.mean - overall.mean
age40to65.diff.se <- sqrt(age40to65.se^2 + overall.se^2)
age65p.diff <- age65p.mean - overall.mean
age65p.diff.se <- sqrt(age65p.se^2 + overall.se^2)

# store estimates 
m <- c(agel40.diff, age40to65.diff, age65p.diff)
s <- c(agel40.diff.se, age40to65.diff.se, age65p.diff.se)
haq.lprog.age <- data.table(est = m, se = s)
haq.lprog.age[, lower := est - qnorm(.975) * se]
haq.lprog.age[, upper := est + qnorm(.975) * se]
rownames(haq.lprog.age) <- c("< 40", "40-64", " >= 65")
save(haq.lprog.age, 
     file = "../data/haq-lprog-age.rda", compress = "bzip2")
