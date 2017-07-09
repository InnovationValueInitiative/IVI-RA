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
               "rituximab", "tofacitinib citrate", "certolizumab pegol", "golimumab")
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
nma.acr.naive.coda[["d[36]"]] <- nma.acr.naive.coda[["d[1]"]]
nma.acr.naive.coda[["d[37]"]] <- NA 
nma.acr.naive.coda[["d[38]"]] <- NA 
nma.acr.naive.crosswalk[, coda_num := paste0("d[", num, "]")]
nma.acr.naive.crosswalk <- nma.acr.naive.crosswalk[sname != ""]
pos <- match(therapy.snames, nma.acr.naive.crosswalk$sname)
pos <- pos[!is.na(pos)]
nma.acr.naive.crosswalk <- nma.acr.naive.crosswalk[pos]
coda.num <- nma.acr.naive.crosswalk$coda_num
nma.acr.naive.coda <- nma.acr.naive.coda[, c("A", "z[1]", "z[2]", "z[3]", coda.num), with = FALSE]
setnames(nma.acr.naive.coda, colnames(nma.acr.naive.coda),
         c("A", "z1", "z2", "z3", paste0("d_", nma.acr.naive.crosswalk$sname)))
nma.acr.naive <- list(mean = apply(nma.acr.naive.coda, 2, mean),
                      vcov = cov(nma.acr.naive.coda))

### bio experienced
# hr assumed to have lower bound of .75 and upper bound of 0.92


# #### NMA FROM ICON
# ### Overall
# ## parameters
# A <- 0.563091333129146
# z2 <- list(mean = 0.62193993520384, sd = 0.00721932018954759)
# z3 <- list(mean = 1.20738138209854, sd = 0.0100302799769332)
# delta <- list(mean =  rep(NA, nther), sd = rep(NA, nther))
# names(delta$mean) <- names(delta$sd) <- therapy.snames 
# delta$mean["cdmards"] <- 0
# delta$sd["cdmards"] <- 0
# delta$mean["abtivmtx"] <- -0.810586730442689
# delta$sd["abtivmtx"] <- 0.102399606960613
# delta$mean["adamtx"] <- -0.773424375324315
# delta$sd["adamtx"] <- 0.0998397767894423
# delta$mean["ada"] <- -0.348079540649955
# delta$sd["ada"] <- 0.233183736623691
# delta$mean["tt"] <- -0.703516704090787
# delta$sd["tt"] <- 0.277000941871632
# delta$mean["etnmtx"] <- -0.87329542983495
# delta$sd["etnmtx"] <- 0.130169676974215
# delta$mean["etn"] <- -0.528787773611509
# delta$sd["etn"] <- 0.151360422595065
# delta$mean["golmtx"] <- -0.777575893107503
# delta$sd["golmtx"] <- 0.12321948919039
# delta$mean["ifxmtx"] <- -0.71150417108933
# delta$sd["ifxmtx"] <- 0.126373150537479
# delta$mean["placebo"] <- 100    # dummy value, set so that ACR < 20 with probability 1 
# delta$sd["placebo"] <- 0
# delta$mean["tczmtx"] <- -0.849976396987633
# delta$sd["tczmtx"] <- 0.242289152342482
# delta$mean["tcz"] <- -0.873409071106643
# delta$sd["tcz"] <- 0.283288795703314
# delta$mean["czpmtx"] <- -1.18034219680414
# delta$sd["czpmtx"] <- 0.110622277385538
# delta$mean["abtscmtx"] <- -0.810586730442689 # NOTE ASSUMED TO BE SAME AS IV!!!!!
# delta$sd["abtscmtx"] <- 0.102399606960613
# delta$mean["nbt"] <- 100 # dummy value, set so that ACR < 20 with probability 1
# delta$sd["nbt"] <- 0
# delta$mean["rtxmtx"] <- -0.794692693812092
# delta$sd["rtxmtx"] <- 0.114214486861968
# delta$mean["tofmtx"] <- -0.754649219587801
# delta$sd["tofmtx"] <- 0.12753876539542
# delta$mean["rtx"] <- -0.602894127529288
# delta$sd["rtx"] <- 0.209130659264918
# delta$mean["tof"] <- -0.147774537169544
# delta$sd["tof"] <- 0.398774886986443
# delta$mean["czp"] <- -0.540154951487203
# delta$sd["czp"] <- 0.37893269604465
# delta$mean["gol"] <- 0.210217681646764
# delta$sd["gol"] <- 0.348248823065528
# 
# ## probability at means
# p <- acr_prob(A, z2$mean, z3$mean, delta$mean)
# 
# ## store results in list
# icon_acr_nma <- list(mean = c(A, z2$mean, z3$mean, delta$mean),
#                      vcov = diag(c(0, z2$sd^2, z3$sd^2, delta$sd^2)),
#                      p = p)
# names(icon_acr_nma$mean) <- c("A", "z2", "z3", names(delta$mean))
# 
# ### TIM naive
# ## parameters
# A <- .5385
# z2 <- list(mean = 0.62171464, sd = 0.0075446477739013)
# z3 <- list(mean = 1.20949933333333, sd = 0.0105001763691831)
# delta <- list(mean =  rep(NA, nther), sd = rep(NA, nther))
# names(delta$mean) <- names(delta$sd) <- therapy.snames 
# delta$mean["cdmards"] <- 0
# delta$sd["cdmards"] <- 0
# delta$mean["abtivmtx"] <- -0.809142111666667
# delta$sd["abtivmtx"] <- 0.123000025214947
# delta$mean["adamtx"] <- -0.769457638333333
# delta$sd["adamtx"] <- 0.1078885028064
# delta$mean["ada"] <- -0.348079540649955
# delta$sd["ada"] <- 0.233183736623691
# delta$mean["tt"] <- -0.703516704090787
# delta$sd["tt"] <- 0.277000941871632
# delta$mean["etnmtx"] <- -0.872833673333333
# delta$sd["etnmtx"] <- 0.121932640740375
# delta$mean["etn"] <- -0.535766297733333
# delta$sd["etn"] <- 0.155748940193373
# delta$mean["golmtx"] <- -0.8133894416666673
# delta$sd["golmtx"] <- 0.147094404471251
# delta$mean["ifxmtx"] <- -0.762230790666667
# delta$sd["ifxmtx"] <- 0.143644467849748
# delta$mean["placebo"] <- 100    # dummy value, set so that ACR < 20 with probability 1 
# delta$sd["placebo"] <- 0
# delta$mean["tczmtx"] <- -0.791568189083333
# delta$sd["tczmtx"] <- 0.265996916175679
# delta$mean["tcz"] <- -0.837710285040667
# delta$sd["tcz"] <- 0.310588145282049
# delta$mean["czpmtx"] <- -1.19416948833333
# delta$sd["czpmtx"] <- 0.119851834884759
# delta$mean["abtscmtx"] <- -0.809142111666667 # NOTE ASSUMED TO BE SAME AS IV!!!!!
# delta$sd["abtscmtx"] <- 0.123000025214947
# delta$mean["nbt"] <- 100 # dummy value, set so that ACR < 20 with probability 1
# delta$sd["nbt"] <- 0
# delta$mean["rtxmtx"] <- -0.680662496166667
# delta$sd["rtxmtx"] <- 0.169679667357375
# delta$mean["tofmtx"] <- -0.748702311666667
# delta$sd["tofmtx"] <- 0.138277342355426
# delta$mean["rtx"] <- -0.602894127529288 # NOTE ASSUMED TO BE SAME AS IN COMBINED CASE
# delta$sd["rtx"] <- 0.209130659264918
# delta$mean["tof"] <- -0.114165030696333
# delta$sd["tof"] <- 0.401271826028596
# delta$mean["czp"] <- -0.488228497576
# delta$sd["czp"] <- 0.396949683345864
# delta$mean["gol"] <- 0.210217681646764 # NOTE ASSUMED TO BE SAME AS IN COMBINED CASE
# delta$sd["gol"] <- 0.348248823065528
# 
# ## probability at means
# p <- acr_prob(A, z2$mean, z3$mean, delta$mean)
# 
# ## store results in list
# icon_acr_nma_naive <- list(mean = c(A, z2$mean, z3$mean, delta$mean),
#                      vcov = diag(c(0, z2$sd^2, z3$sd^2, delta$sd^2)),
#                      p = p)
# names(icon_acr_nma_naive$mean) <- c("A", "z2", "z3", names(delta$mean))
# 
# ### TIM experienced
# ## parameters
# A <- .7454
# z2 <- list(mean = 0.616629345, sd = 0.025710244721687)
# z3 <- list(mean = 1.11857970333333, sd = 0.036565591741695)
# delta <- list(mean =  icon_acr_nma_naive$mean[-c(1:3)], 
#               sd = sqrt(diag(icon_acr_nma_naive$vcov))[-c(1:3)])
# names(delta$sd) <- names(delta$mean)
# delta$mean["abtivmtx"] <- 0.68978946169
# delta$sd["abtivmtx"] <- 2.28054078395492
# delta$mean["abtscmtx"] <- 0.68978946169
# delta$sd["abtscmtx"] <- 2.28054078395492
# delta$mean["golmtx"] <- -1.238279128404
# delta$sd["golmtx"] <- 2.69756416402609
# delta$mean["rtxmtx"] <- -1.16267591078228
# delta$sd["rtxmtx"] <- 2.98659459316976
# 
# ## probability at means
# p <- acr_prob(A, z2$mean, z3$mean, delta$mean)
# 
# ## store results in list
# icon_acr_nma_exp <- list(mean = c(A, z2$mean, z3$mean, delta$mean),
#                            vcov = diag(c(0, z2$sd^2, z3$sd^2, delta$sd^2)),
#                            p = p)
# names(icon_acr_nma_exp$mean) <- c("A", "z2", "z3", names(delta$mean))

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
cost <- data.table(read_excel("cost.xlsx", sheet = "RData"))
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
therapy.pars <- list(info = therapy.info,
                    nma.acr.naive = nma.acr.naive,
                    nice.nma.acr.naive = nice.nma.acr.naive,
                    haq.lprog = haq.lprog, cost = cost,
                    si = si)
save(therapy.pars, file = "../data/therapy-pars.rda", compress = "bzip2")
