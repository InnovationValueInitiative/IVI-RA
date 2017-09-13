rm(list = ls())
library("data.table")

# TREATMENT ATTRIBUTES --------------------------------------------------------
load("../data/treatments.rda")
utility.tx.attr <- list()

# design matrix
utility.tx.attr$x <- model.matrix(~ route + years_since_approval,treatments)[, -1]
rownames(utility.tx.attr$x) <- treatments$sname

# parameters
utility.tx.attr$coef <- data.table(var = colnames(utility.tx.attr$x),
                                   lower = rep(0, ncol(utility.tx.attr$x)),
                                   upper = rep(0, ncol(utility.tx.attr$x)))
# save
save(utility.tx.attr, file = "../data/utility-tx-attr.rda", compress = "bzip2") 

# ALAVA (2013) MIXTURE MODEL --------------------------------------------------
nice <- read.csv("nice-painhaq.csv")

## pain
# parameters from Sarzi-Puttini et al (2002)
pain.mean <- 61.65
haq.mean <- 1.39
interval <- 19.10
pain.var <- 19.10^2
haq.var <- .59^2
painhaq.cor <- .421

# example sampling
pain_sample <- function(nsims, haq){
  pain.condmean <- pain.mean + painhaq.cor * sqrt((pain.var/haq.var)) * (haq - haq.mean)
  pain.condsd <- sqrt((1 - painhaq.cor^2) * pain.var)
  pain.sample <- rnorm(nsims, pain.condmean, pain.condsd)
  return(pain.sample)
}
tmp <- pain_sample(1000, 1.5)
hist(tmp)
summary(tmp)

# correlation implied by nice
nice.lm <- lm(pain ~ haq, nice)
nice.cor <- coef(nice.lm)[2] * sqrt(haq.var/pain.var)
nice.lm2 <- lm(pain ~ haq + I(haq^2) , nice)

# save
pain <- list(pain.mean = pain.mean, haq.mean = haq.mean, 
             pain.var = pain.var, haq.var = haq.var,
             painhaq.cor = nice.cor)
save(pain, 
     file = "../data/pain.rda", compress = "bzip2")

# utility
coef <- read.csv("coef-alava-2013.csv")
vcov <- read.csv("vcov-alava-2013.csv")
vcov.names <- vcov[, 1]
vcov <- as.matrix(vcov[, -1])
rownames(vcov) <- vcov.names
# index <- list()
# delta.names <- paste0("delta", c(1, 2, 3))
# parameters <- unique(coef$parameter)
# parameters <- parameters[!parameters %in% delta.names]
# for (par in parameters){
#   index[[par]] <- which(coef$parameter == par)
# }
# index[["delta"]] <- which(coef$parameter %in% delta.names)
utility.mixture <- list(coef = coef, vcov = vcov)
save(utility.mixture, 
     file = "../data/utility-mixture.rda", compress = "bzip2")

# WAILOO (2006) LOGIT MODEL ---------------------------------------------------
est <- c(2.0734, .0058, .0023, -.2004, -.2914, .0249, -.8647)
se <- c(.0263, .0004, .0004, .0101, .0118, .0028, .0103)
vars <- c("int", "age", "dis_dur", "haq0", "male", "prev_dmards", "haq")
utility.wailoo <- data.table(var = vars, est = est, se = se)
save(utility.wailoo, 
     file = "../data/utility-wailoo.rda", compress = "bzip2")