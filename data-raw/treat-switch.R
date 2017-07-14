rm(list = ls())
library("data.table")
source("func.R")

# convert odds ratio to log odds ratio
or <- c(1.94, 3.39)
or.l <- c(1.36, 2.27)
or.u <- c(2.75, 5.06)
logor <- log(or) 
logor.l <- log(or.l)
logor.u <- log(or.u)
logor.se <- c(se_normal(logor.l[1], logor.u[1]), se_normal(logor.l[2], logor.u[2]))

# probability by disease activity level
p.moderate <- .163 # from CORRONA
logor.int <- logit(p.moderate) - logor[1]
or.int <- exp(logor.int)
logor.int.l <- logor.int - .03
logor.int.u <- logor.int + .03
or.int.l <- exp(logor.int.l)
or.int.u <- exp(logor.int.u)
logor.int.se <- se_normal(logor.int.l, logor.int.u)
p.low <- 1/(1 + exp(-logor.int))
p.high <- 1/(1 + exp(-logor.int - logor[2]))

# beta distribution for probabilities for four categories of disease activity
# remission, low, moderate, high
n.low <- 44 + 293
n.moderate <- 76 + 551
n.high <- 110 + 679
p <- c(p.low, p.low, p.moderate, p.high)
n <- c(round(n.low/2), round(n.low/2), n.moderate, n.high)
alpha <- round(p * n)
beta <- round(n - alpha)
p.mean <- alpha/(alpha + beta)
p.se <- sqrt((alpha * beta)/((alpha + beta)^2 * (alpha + beta + 1)))
p.lower <- qbeta(.025, alpha, beta)
p.upper <- qbeta(.975, alpha, beta)

# save
treat.switch <- data.table(disease_activity = c("Remission", "Low", "Moderate", "High"),
                           p = p.mean, p_se = p.se,
                           p_lower = p.lower, p_upper = p.upper,
                           alpha = alpha, beta = beta,
                           n = n)
save(treat.switch, file = "../data/treat-switch.rda", compress = "bzip2")
