rm(list = ls())
library("data.table")
source("func.R")

# convet odds ratio to log odds ratio
or <- c(1.94, 3.39)
or.l <- c(1.36, 2.27)
or.u <- c(2.75, 5.06)
logor <- log(or) 
logor.l <- log(or.l)
logor.u <- log(or.u)
logor.se <- c(se_normal(logor.l[1], logor.u[1]), se_normal(logor.l[2], logor.u[2]))

# baseline probability from CORRONA
p.moderate <- .163
logor.int <- logit(p.moderate) - logor[1]
or.int <- exp(logor.int)
logor.int.l <- logor.int - .03
logor.int.u <- logor.int + .03
or.int.l <- exp(logor.int.l)
or.int.u <- exp(logor.int.u)
logor.int.se <- se_normal(logor.int.l, logor.int.u)
p.low <- 1/(1 + exp(-logor.int))


# save
treat.switch <- data.table(variable = c("intercept", "da_moderate", "da_high"), 
                           or = c(or.int, or), 
                           or_lower = c(or.int.l, or.l), or_upper = c(or.int.u, or.u),
                           logor = c(logor.int, logor), logor_se = c(logor.int.se, logor.se), 
                           logor_lower = c(logor.int.l, logor.l), 
                           logor_upper = c(logor.int.u, logor.u))
save(treat.switch, file = "../data/treat-switch.rda", compress = "bzip2")
