rm(list = ls())
library("data.table")
acr.cats <- c("ACR < 20", "ACR 20-50", "ACR 50-70", "ACR 70+")

# ACR TO EULAR -----------------------------------------------------------------
acr2eular <- matrix(c(755, 4, 2, 0, 136, 27, 2, 2, 57, 26, 10, 2), 
                    nrow = 4, ncol = 3, byrow = FALSE)
rownames(acr2eular) <- acr.cats
colnames(acr2eular) <- c("eular_none", "eular_moderate", "eular_good")
save(acr2eular, file = "../data/acr2eular.rda", compress = "bzip2")

# ACR TO HAQ -------------------------------------------------------------------
acr2haq <- data.table(acr = acr.cats,
                      mean = c(-.11, -.44, -.76, -1.07),
                      se = c(.06765, .05657, .09059, .07489))
save(acr2haq, file = "../data/acr2haq.rda", compress = "bzip2")

# ACR TO SDAI/CDAI/DAS28 -------------------------------------------------------
## calculate mean in group 1 given overall mean, n in each group, and mean in group 2
grp_mean <- function(n1, n2, mean2, mean_all){
  n <- n1 + n2
  p1 <- n1/n; p2 <- n2/n
  mean1 <- (mean_all - p2 * mean2)/p1
  # mean_all <- p1 * mean1 + p2 * mean2 
  return(mean1)
}

## leflunomide dataset and inception cohort
acr.lef.n <- c(948, 398, 342, 151)
acr.inc.n <- c(28, 9, 11, 43)

## acr to sdai
# leflunomide dataset 
acr70.sdai.lef <- -41
acr50.sdai.lef <- grp_mean(acr.lef.n[3], acr.lef.n[4], -41, -37) 
acr20.sdai.lef <- grp_mean(acr.lef.n[2], acr.lef.n[3] + acr.lef.n[4], -37, -34) 
stopifnot(acr.lef.n[2]/sum(acr.lef.n[-1]) * acr20.sdai.lef + 
            acr.lef.n[3]/sum(acr.lef.n[-1]) * acr50.sdai.lef +
            acr.lef.n[4]/sum(acr.lef.n[-1]) * acr70.sdai.lef == -34)

# inception cohort
acr70.sdai.inc <- -30.1
acr50.sdai.inc <- grp_mean(acr.inc.n[3], acr.inc.n[4], -30.1, -27) 
acr20.sdai.inc <- grp_mean(acr.inc.n[2], acr.inc.n[3] + acr.inc.n[4], -27, -25.1) 
stopifnot(round(acr.inc.n[2]/sum(acr.inc.n[-1]) * acr20.sdai.inc + 
            acr.inc.n[3]/sum(acr.inc.n[-1]) * acr50.sdai.inc +
            acr.inc.n[4]/sum(acr.inc.n[-1]) * acr70.sdai.inc, 1) == -25.1)

## acr to cdai
# inception cohort
acr70.cdai.inc <- -27.6
acr50.cdai.inc <- grp_mean(acr.inc.n[3], acr.inc.n[4], -27.6, -24.6) 
acr20.cdai.inc <- grp_mean(acr.inc.n[2], acr.inc.n[3] + acr.inc.n[4], -24.6, -22.7) 
stopifnot(round(acr.inc.n[2]/sum(acr.inc.n[-1]) * acr20.cdai.inc + 
            acr.inc.n[3]/sum(acr.inc.n[-1]) * acr50.cdai.inc +
            acr.inc.n[4]/sum(acr.inc.n[-1]) * acr70.cdai.inc, 1) == -22.7)

## acr to das28
# inception cohort
acr70.das28.inc <- -3.31
acr50.das28.inc <- grp_mean(acr.inc.n[3], acr.inc.n[4], -3.31, -2.95) 
acr20.das28.inc <- grp_mean(acr.inc.n[2], acr.inc.n[3] + acr.inc.n[4], -2.95, -2.75) 
stopifnot(round(acr.inc.n[2]/sum(acr.inc.n[-1]) * acr20.das28.inc + 
                  acr.inc.n[3]/sum(acr.inc.n[-1]) * acr50.das28.inc +
                  acr.inc.n[4]/sum(acr.inc.n[-1]) * acr70.das28.inc, 2) == -2.75)


## save
acr2sdai <- list(leflunomide = data.table(acr = acr.cats,
                                        mean = c(0, acr20.sdai.lef,
                                                 acr50.sdai.lef, acr70.sdai.lef)),
                 inception = data.table(acr = acr.cats,
                                        mean = c(0, acr20.sdai.inc,
                                                 acr50.sdai.inc, acr70.sdai.inc)))
acr2sdai$leflunomide[, lower := mean * 1.2]
acr2sdai$leflunomide[, upper := mean * .8]
acr2sdai$inception[, lower := mean * 1.2]
acr2sdai$inception[, upper := mean * .8]
acr2cdai <- list(inception = data.table(acr = acr.cats,
                                        mean = c(0, acr20.cdai.inc,
                                                 acr50.cdai.inc, acr70.cdai.inc)))
acr2cdai$inception[, lower := mean * 1.2]
acr2cdai$inception[, upper := mean * .8]
acr2das28 <- list(inception = data.table(acr = acr.cats,
                                        mean = c(0, acr20.das28.inc,
                                                 acr50.das28.inc, acr70.das28.inc)))
acr2das28$inception[, lower := mean * 1.2]
acr2das28$inception[, upper := mean * .8]
save(acr2sdai, file = "../data/acr2sdai.rda", compress = "bzip2")
save(acr2cdai, file = "../data/acr2cdai.rda", compress = "bzip2")
save(acr2das28, file = "../data/acr2das28.rda", compress = "bzip2")