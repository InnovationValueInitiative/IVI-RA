rm(list = ls())
library("iviRA")
library("data.table")
source("func.R")
treatments <- fread("treatments.csv")

# HAQ PROGRESSION RATE BY TREATMENT --------------------------------------------
# progression rate with and without biologics from Wolfe & Michaud (2010)
hlp.nb <- .031
hlp.nb.se <- se_normal(.026, .036)
hlp.b <- -.001
hlp.b.se <- se_normal(-.004, .002)

# progression rate for biologic specific from Michaud (2011)
hlp.ifxmtx <- -.001
hlp.ifxmtx.se <- se_normal(-.008, .006) 
hlp.etn <- -.003
hlp.etn.se <- se_normal(-.010, .005)  
hlp.etnmtx <- -.005
hlp.etnmtx.se <- se_normal(-.012, .001)  
hlp.adamtx <- -.003
hlp.adamtx.se <- se_normal(-.018, .011)  
hlp.ada <- .012
hlp.ada.se <- se_normal(-.002, .025)  

# biologic specific progression rates from Michaud (2011) and above
haq.lp.tx <- data.table(sname = treatments$sname, est = 0, se = 0)
haq.lp.tx[sname %in% c("cdmards", "placebo","nbt"), 
         c("est", "se") :=  list(hlp.nb, hlp.nb.se)]
haq.lp.tx[sname %in% c("abtivmtx", "golmtx", "tcz", "tczmtx",
                        "czpmtx", "abtscmtx","rtxmtx", "tofmtx", 
                      "rtx", "tof", "czp", "gol"), 
         c("est", "se") :=  list(hlp.b, hlp.b.se)]
haq.lp.tx[sname == "ifxmtx", 
         c("est", "se") := list(hlp.ifxmtx, hlp.ifxmtx.se)]
haq.lp.tx[sname == "etn", 
         c("est", "se") := list(hlp.etn, hlp.etn.se)]
haq.lp.tx[sname == "etnmtx", 
         c("est", "se") := list(hlp.etnmtx, hlp.etnmtx.se)]
haq.lp.tx[sname == "ada", 
         c("est", "se") := list(hlp.ada, hlp.ada.se)]
haq.lp.tx[sname == "adamtx", 
         c("est", "se") := list(hlp.adamtx, hlp.adamtx.se)]

# confidence bounds
haq.lp.tx[, lower := est - qnorm(.975) * se]
haq.lp.tx[, upper := est + qnorm(.975) * se]

# HAQ PROGRESSION RATE BY AGE --------------------------------------------------
# results from table 1 in Michaud (2011)
haq.lp.overall.mean <- .014
haq.lp.overall.lower <- .012
haq.lp.overall.upper <- .015

haq.lp.agel40.mean <- -.006
haq.lp.agel40.lower <- -.014
haq.lp.agel40.upper <- .002

haq.lp.age40to65.mean <- .006
haq.lp.age40to65.lower <- .004
haq.lp.age40to65.upper <- .007

haq.lp.age65p.mean <- .031
haq.lp.age65p.lower <- .011
haq.lp.age65p.upper <- .019

# calculate standard errors 
haq.lp.overall.se <- se_normal(haq.lp.overall.lower, haq.lp.overall.upper)
haq.lp.agel40.se <- se_normal(haq.lp.agel40.lower, haq.lp.agel40.upper)
haq.lp.age40to65.se <- se_normal(haq.lp.age40to65.lower, haq.lp.age40to65.upper)
haq.lp.age65p.se <- se_normal(haq.lp.age65p.lower, haq.lp.age65p.upper)

# differences
dhaq.lp.agel40 <- haq.lp.agel40.mean - haq.lp.overall.mean
dhaq.lp.agel40.se <- sqrt(haq.lp.agel40.se^2 + haq.lp.overall.se^2)
dhaq.lp.age40to65 <- haq.lp.age40to65.mean - haq.lp.overall.mean
dhaq.lp.age40to65.se <- sqrt(haq.lp.age40to65.se^2 + haq.lp.overall.se^2)
dhaq.lp.age65p <- haq.lp.age65p.mean - haq.lp.overall.mean
dhaq.lp.age65p.se <- sqrt(haq.lp.age65p.se^2 + haq.lp.overall.se^2)

# store estimates 
dhaq.lp.age.mean <- c(dhaq.lp.agel40, dhaq.lp.age40to65, dhaq.lp.age65p)
dhaq.lp.age.se <- c(dhaq.lp.agel40.se, dhaq.lp.age40to65.se, dhaq.lp.age65p.se)
dhaq.lp.age <- data.table(age = c("< 40", "40-64", " >= 65"),
                         est = dhaq.lp.age.mean, se = dhaq.lp.age.se)
dhaq.lp.age[, lower := est - qnorm(.975) * se]
dhaq.lp.age[, upper := est + qnorm(.975) * se]

# SAVE PARAMETERS --------------------------------------------------------------
haq.lprog <- list(tx = haq.lp.tx, diff.age = dhaq.lp.age)
save(haq.lprog, file = "../data/haq-lprog.rda", compress = "bzip2")
