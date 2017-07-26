rm(list = ls())
library("iviRA")
library("data.table")
library("MCMCpack")
source("func.R")
treatments <- fread("treatments.csv")
nther <- nrow(treatments)

# LOAD DATA --------------------------------------------------------------------
nma.acr.naive.coda <- fread("nma-acr-naive-re-coda.csv")
nma.acr.naive.crosswalk <- fread("nma-acr-naive-re-crosswalk.csv")
nma.haq.naive.coda <- fread("nma-haq-naive-re-coda.csv")
nma.haq.naive.crosswalk <- fread("nma-haq-naive-re-crosswalk.csv")
nma.das28.naive.coda <- fread("nma-das28-naive-re-coda.csv")
nma.das28.naive.crosswalk <- fread("nma-das28-naive-re-crosswalk.csv")

# ACR RESPONSE PROBABILITIES ---------------------------------------------------
#### NMA FROM NICE
### bio naive
## table 37 Stevenson (also assuming RTX is equal to 
# i.v. ABT per p. 248 in NICE report)
po <- vector(mode = "list", nther)
names(po) <- treatments$sname 
po[["cdmards"]] <- c(0.298, 0.123, 0.042)
po[["abtivmtx"]] <- c(0.573, 0.328, 0.156)
po[["adamtx"]] <- c(0.615, 0.368, 0.183)
po[["ada"]] <- c(0.499, 0.264, 0.115)
po[["tt"]] <- c(0.503, 0.266, 0.117)
po[["etnmtx"]] <- c(0.713, 0.472, 0.263)
po[["etn"]] <- c(0.645, 0.398, 0.205)
po[["golmtx"]] <- c(0.642, 0.395, 0.202)
po[["ifxmtx"]] <- c(0.595, 0.348, 0.169)
po[["placebo"]] <- c(0.175, 0.059, 0.016)
po[["tczmtx"]] <- c(0.706, 0.464, 0.256)
po[["tcz"]] <- c(0.717, 0.477, 0.266)
po[["czpmtx"]] <- c(0.564, 0.319, 0.150)
po[["abtscmtx"]] <- c(0.638, 0.391, 0.199)
po[["nbt"]] <- c(0, 0, 0)
po[["rtxmtx"]] <- c(0.573, 0.328, 0.156) # assumed to be same as abt iv
po[["tofmtx"]] <- c(NA, NA, NA)
po[["rtx"]] <- c(NA, NA, NA)
po[["tof"]] <- c(NA, NA, NA)
po[["czp"]] <- c(NA, NA, NA)
po[["gol"]] <- c(NA, NA, NA)
po <- do.call("rbind", po)
po <- cbind(1 - po[, 1], po)

## probabilities in mutually exclusive categories
p <- matrix(NA, nrow = nther, ncol = 4)
rownames(p) <- rownames(po)
p[, 1] <- 1 - po[, 2] 
p[, 2] <- po[, 2] - po[, 3]
p[, 3] <- po[, 3] - po[, 4]
p[, 4] <- po[, 4]

## Guess Parameters of Dirichlet Distribution
# simulate cDMARDs
nsims <- 1000
nobs <- 500
acr.probs <- MCMCpack::rdirichlet(1000, nobs * p[which(rownames(p) == "cdmards"), ])
acr.probs2 <- matrix(NA, nrow = nsims, ncol = 3)
acr.probs2[, 1] <- 1 - acr.probs[, 1]
acr.probs2[, 2] <- acr.probs[, 3] + acr.probs[, 4]
acr.probs2[, 3] <- acr.probs[, 4]

# check 95% credible interval for cDMARds 
# In Steveson, ACR 20 = .298 (.255-.344), ACR 50 = .123 (.098-.1530),
# ACR 70 = .042 (.031-.056)
apply(acr.probs2, 2, quantile, c(.025, .975))

## Store Results
nma.acr.naive.nice <- list(p = p, p.overlap = po, nobs = nobs)

#### IVI NMA
### bio naive
## missing trials
nma.acr.naive.coda <- nma_add_missing(nma.acr.naive.coda, nma.acr.naive.crosswalk)
nma.acr.naive.crosswalk[, coda_num := paste0("d[", num, "]")]
nma.acr.naive.crosswalk <- nma.acr.naive.crosswalk[sname != ""]

# match nma therapies with model therapies
nma.acr.naive.crosswalk <- dt_reorder(nma.acr.naive.crosswalk, "sname",
                                      treatments$sname)

# reorder coda 
coda.num <- nma.acr.naive.crosswalk$coda_num
nma.acr.naive.coda <- nma.acr.naive.coda[, c("A", "z[1]", "z[2]", "z[3]", coda.num), with = FALSE]
setnames(nma.acr.naive.coda, colnames(nma.acr.naive.coda),
         c("A", "z1", "z2", "z3", paste0("d_", nma.acr.naive.crosswalk$sname)))
nma.acr.naive <- list(mean = apply(nma.acr.naive.coda, 2, mean),
                      vcov = cov(nma.acr.naive.coda))

# CHANGE IN HAQ ----------------------------------------------------------------
#### IVI NMA
### bio naive
## missing trials
nma.haq.naive.coda <- nma_add_missing(nma.haq.naive.coda, nma.haq.naive.crosswalk)
nma.haq.naive.crosswalk[, coda_num := paste0("d[", num, "]")]
nma.haq.naive.crosswalk <- nma.haq.naive.crosswalk[sname != ""]

# match nma therapies with model therapies
nma.haq.naive.crosswalk <- dt_reorder(nma.haq.naive.crosswalk, "sname",
                                      treatments$sname)

# reorder coda
coda.num <- nma.haq.naive.crosswalk$coda_num
nma.haq.naive.coda <- nma.haq.naive.coda[, c("A", coda.num), with = FALSE]
setnames(nma.haq.naive.coda, colnames(nma.haq.naive.coda),
         c("A", paste0("d_", nma.haq.naive.crosswalk$sname)))
nma.haq.naive <- list(mean = apply(nma.haq.naive.coda, 2, mean),
                      vcov = cov(nma.haq.naive.coda))

# CHANGE IN DAS28 --------------------------------------------------------------
### bio naive
## missing trials
nma.das28.naive.coda <- nma_add_missing(nma.das28.naive.coda, nma.das28.naive.crosswalk)
nma.das28.naive.crosswalk[, coda_num := paste0("d[", num, "]")]
nma.das28.naive.crosswalk <- nma.das28.naive.crosswalk[sname != ""]

# match nma therapies with model therapies
nma.das28.naive.crosswalk <- dt_reorder(nma.das28.naive.crosswalk, "sname",
                                      treatments$sname)

# reorder coda
coda.num <- nma.das28.naive.crosswalk$coda_num
nma.das28.naive.coda <- nma.das28.naive.coda[, c("A", coda.num), with = FALSE]
setnames(nma.das28.naive.coda, colnames(nma.das28.naive.coda),
         c("A", paste0("d_", nma.das28.naive.crosswalk$sname)))
nma.das28.naive <- list(mean = apply(nma.das28.naive.coda, 2, mean),
                      vcov = cov(nma.das28.naive.coda))

# SAVE PARAMETERS --------------------------------------------------------------
# We currently don't have NMA results for triple therapy
tt.indx <- which(treatments$sname == "tt")
nma.acr.naive$mean <- nma.acr.naive$mean[-c(tt.indx + 4)]
nma.acr.naive$vcov <- nma.acr.naive$vcov[-c(tt.indx + 4), -c(tt.indx + 4)]
nma.haq.naive$mean <- nma.haq.naive$mean[-c(tt.indx + 1)]
nma.haq.naive$vcov <- nma.haq.naive$vcov[-c(tt.indx + 1), -c(tt.indx + 1)]
nma.das28.naive$mean <- nma.das28.naive$mean[-c(tt.indx + 1)]
nma.das28.naive$vcov <- nma.das28.naive$vcov[-c(tt.indx + 1), -c(tt.indx + 1)]
nma.acr.naive.nice$p <- nma.acr.naive.nice$p[-tt.indx, ]
nma.acr.naive.nice$p.overlap <- nma.acr.naive.nice$p.overlap[-tt.indx, ]

# save
save(nma.acr.naive, file = "../data/nma-acr-naive.rda", compress = "bzip2")
save(nma.acr.naive.nice, file = "../data/nma-acr-naive-nice.rda", compress = "bzip2")
save(nma.das28.naive, file = "../data/nma-das28-naive.rda", compress = "bzip2")
save(nma.haq.naive, file = "../data/nma-haq-naive.rda", compress = "bzip2")
