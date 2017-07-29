rm(list = ls())
library("iviRA")
library("MCMCpack")
library("readxl")
library("data.table")
source("func.R")
therapy.names <- c("cDMARDs", "abatacept IV + methotrexate", "adalimumab + methotrexate",
               "adalimumab", "Triple therapy", "etanercept + methotrexate", "etanercept", 
               "golimumab + methotrexate", "infliximab + methotrexate", "placebo",
               "tocilizumab + methotrexate", "tocilizumab", "certolizumab pegol + methotrexate",
               "abatacept SC + methotrexate", "non-biologic", 
               "rituximab + methotrexate", "tofacitinib citrate + methotrexate",
               "rituximab", "tofac itinib citrate", "certolizumab pegol", "golimumab")
therapy.mnames <- c("cDMARDs", "ABT IV + MTX", "ADA + MTX", "ADA", "Triple therapy",
                    "ETN + MTX", "ETN", "GOL + MTX", "IFX + MTX", "Placebo", "TCZ + MTX",
                    "TCZ", "CZP + MTX", "ABT SC + MTX", "NBT", "RTX + MTX", "TOF + MTX",
                    "RTX", "TOF", "CZP", "GOL")
therapy.snames <- c("cdmards", "abtivmtx", "adamtx", "ada", "tt", "etnmtx",
                    "etn", "golmtx", "ifxmtx", "placebo", "tczmtx", "tcz", 
                    "czpmtx", "abtscmtx", "nbt", "rtxmtx", "tofmtx", "rtx", "tof", 
                    "czp", "gol")
therapy.info <- data.table(name = therapy.names, mname = therapy.mnames, 
                           sname = therapy.snames)
therapy.info[, biologic := ifelse(sname %in% c("cdmards", "tt", "placebo",
                                               "nbt"), 0, 1)]

nther <- nrow(therapy.info)

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
po <- vector(mode = "list", length(therapy.snames))
names(po) <- therapy.snames 
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
acr.probs <- rdirichlet(1000, nobs * p[which(rownames(p) == "cdmards"), ])
acr.probs2 <- matrix(NA, nrow = nsims, ncol = 3)
acr.probs2[, 1] <- 1 - acr.probs[, 1]
acr.probs2[, 2] <- acr.probs[, 3] + acr.probs[, 4]
acr.probs2[, 3] <- acr.probs[, 4]

# check 95% credible interval for cDMARds 
# In Steveson, ACR 20 = .298 (.255-.344), ACR 50 = .123 (.098-.1530),
# ACR 70 = .042 (.031-.056)
apply(acr.probs2, 2, quantile, c(.025, .975))

## Store Results
nice.nma.acr.naive <- list(p = p, p.overlap = po, nobs = nobs)

#### IVI NMA
### bio naive
## missing trials
nma.acr.naive.coda <- nma_add_missing(nma.acr.naive.coda, nma.acr.naive.crosswalk)
nma.acr.naive.crosswalk[, coda_num := paste0("d[", num, "]")]
nma.acr.naive.crosswalk <- nma.acr.naive.crosswalk[sname != ""]

# match nma therapies with model therapies
nma.acr.naive.crosswalk <- dt_reorder(nma.acr.naive.crosswalk, "sname",
                                      therapy.snames)

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
                                      therapy.snames)

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
                                      therapy.snames)

# reorder coda
coda.num <- nma.das28.naive.crosswalk$coda_num
nma.das28.naive.coda <- nma.das28.naive.coda[, c("A", coda.num), with = FALSE]
setnames(nma.das28.naive.coda, colnames(nma.das28.naive.coda),
         c("A", paste0("d_", nma.das28.naive.crosswalk$sname)))
nma.das28.naive <- list(mean = apply(nma.das28.naive.coda, 2, mean),
                      vcov = cov(nma.das28.naive.coda))

# HAQ PROGRESSION RATE ---------------------------------------------------------
# progression rate with and without biologics from Wolfe & Michaud (2010)
hp.nb <- .031
hp.nb.se <- se_normal(.026, .036)
hp.b <- -.001
hp.b.se <- se_normal(-.004, .002)

# progression rate for biologic specific from Michaud (2011)
hp.ifxmtx <- -.001
hp.ifxmtx.se <- se_normal(-.008, .006) 
hp.etn <- -.003
hp.etn.se <- se_normal(-.010, .005)  
hp.etnmtx <- -.005
hp.etnmtx.se <- se_normal(-.012, .001)  
hp.adamtx <- -.003
hp.adamtx.se <- se_normal(-.018, .011)  
hp.ada <- .012
hp.ada.se <- se_normal(-.002, .025)  

# biologic specific progression rates from Michaud (2011) and above
haq.lprog <- data.table(sname = therapy.snames, est = 0, se = 0)
haq.lprog[sname %in% c("cdmards", "tt", "placebo","nbt"), 
         c("est", "se") :=  list(hp.nb, hp.nb.se)]
haq.lprog[sname %in% c("abtivmtx", "golmtx", "tcz", "tczmtx",
                        "czpmtx", "abtscmtx","rtxmtx", "tofmtx", 
                      "rtx", "tof", "czp", "gol"), 
         c("est", "se") :=  list(hp.b, hp.b.se)]
haq.lprog[sname == "ifxmtx", 
         c("est", "se") := list(hp.ifxmtx, hp.ifxmtx.se)]
haq.lprog[sname == "etn", 
         c("est", "se") := list(hp.etn, hp.etn.se)]
haq.lprog[sname == "etnmtx", 
         c("est", "se") := list(hp.etnmtx, hp.etnmtx.se)]
haq.lprog[sname == "ada", 
         c("est", "se") := list(hp.ada, hp.ada.se)]
haq.lprog[sname == "adamtx", 
         c("est", "se") := list(hp.adamtx, hp.adamtx.se)]

# confidence bounds
haq.lprog[, lower := est - qnorm(.975) * se]
haq.lprog[, upper := est + qnorm(.975) * se]

# TREATMENT COSTS --------------------------------------------------------------
cost <- data.table(read_excel("cost.xlsx", sheet = "Cost"))
cost[, dose_escalation := .025 * .188]
cost[sname == "ada", dose_escalation := .096 * .912]
cost[sname == "ifx", dose_escalation := .35 * .295]

# ADVERSE EVENTS ---------------------------------------------------------------
### Parameters from Cochrane review (Singh 2011)
si <- list()
si$exp$est <- ifelse(therapy.snames %in% c("cdmards", "tt"), log(.026),
                      ifelse(therapy.snames %in% c("placebo", "nbt"), 
                             log(0 + 10e-11), log(.035)))
names(si$exp$est) <- therapy.snames
si.exp.lrate.se <- se_normal(log(.027), log(.046), .975)
si$exp$vcov <- diag(si.exp.lrate.se^2, length(si$exp$est))
rownames(si$exp$vcov) <- colnames(si$exp$vcov) <- therapy.snames
si$exp$loc.index <- seq_len(length(si$exp$est))
si$exp$anc1.index <- si$exp$anc2.index <- NA

# SAVE PARAMETERS --------------------------------------------------------------
# We currently don't have NMA results for triple therapy
tt.indx <- which(therapy.snames == "tt")
therapy.info <- therapy.info[sname != "tt"]
nma.acr.naive$mean <- nma.acr.naive$mean[-c(tt.indx + 4)]
nma.acr.naive$vcov <- nma.acr.naive$vcov[-c(tt.indx + 4), -c(tt.indx + 4)]
nma.haq.naive$mean <- nma.haq.naive$mean[-c(tt.indx + 1)]
nma.haq.naive$vcov <- nma.haq.naive$vcov[-c(tt.indx + 1), -c(tt.indx + 1)]
nma.das28.naive$mean <- nma.das28.naive$mean[-c(tt.indx + 1)]
nma.das28.naive$vcov <- nma.das28.naive$vcov[-c(tt.indx + 1), -c(tt.indx + 1)]
nice.nma.acr.naive$p <- nice.nma.acr.naive$p[-tt.indx, ]
nice.nma.acr.naive$p.overlap <- nice.nma.acr.naive$p.overlap[-tt.indx, ]
haq.lprog <- haq.lprog[-tt.indx]
si$exp$est <- si$exp$est[-tt.indx]
si$exp$vcov <- si$exp$vcov[-tt.indx, -tt.indx]
si$exp$loc.index <- seq(1, length(si$exp$est))

# save
therapy.pars <- list(info = therapy.info,
                    nma.acr.naive = nma.acr.naive,
                    nma.haq.naive = nma.haq.naive,
                    nma.das28.naive = nma.das28.naive,
                    nice.nma.acr.naive = nice.nma.acr.naive,
                    haq.lprog = haq.lprog, cost = cost,
                    si = si)
save(therapy.pars, file = "../data/therapy-pars.rda", compress = "bzip2")
