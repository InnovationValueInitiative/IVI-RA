rm(list = ls())
library("survival")
library("flexsurv")
library("ggplot2")
library("MASS")
library("data.table")
source("func.R")
theme_set(theme_bw())

# DISEASE ACTIVITY ODDS RATIOS FROM ZHANG (2011) -------------------------------
da.or <- c(1, 1.94, 3.39)

# readjust odds ratio so moderate is the index category
da.or[1] <- da.or[1]/da.or[2] 
da.or[3] <- da.or[3]/da.or[2]
da.or[2] <- 1

# LOAD DIGITIZED SURVIVAL CURVES AND NUMBERS AT RISK ---------------------------
surv.b.g <- read.csv("ttd-bsrbr-eular-good-ggamma.csv")
risk.b.g <- read.csv("ttd-bsrbr-eular-good-risk.csv")
surv.b.m <- read.csv("ttd-bsrbr-eular-moderate-ggamma.csv")
risk.b.m <- read.csv("ttd-bsrbr-eular-moderate-risk.csv")
surv.c <- read.csv("ttd-corrona-allpat.csv")
risk.c <- read.csv("ttd-corrona-allpat-risk.csv")

# CONVERT CURVES TO INDIVIDUAL PATIENT LEVEL PATIENT DATA ----------------------
# algorithm
ipd.b.g <- reconstruct_ipd(surv.b.g, risk.b.g, 1)
ipd.b.g$t <- (ipd.b.g$t/30.42) * 1000 # convert time from 1000s of days to months
ipd.b.g$eular <- "Good"
ipd.b.m <- reconstruct_ipd(surv.b.m, risk.b.m, 1)
ipd.b.m$eular <- "Moderate"
ipd.b.m$t <- (ipd.b.m$t/30.42) * 1000
ipd.c <- reconstruct_ipd(surv.c, risk.c, 1, tot_events = 3584)
ipd.b <- rbind(ipd.b.m, ipd.b.g)

# check kaplan meier of curves
## moderate responders
km.b.m <- survfit(Surv(t, event) ~1,
                       data = ipd.b.m)
p.km.b.m <- km_plot(km.b.m)
ggsave("figs/ttd-bsrbr-ipdsurv-eular-moderate-km.pdf",
       p.km.b.m, height = 5, width = 7)

## good responders
km.b.g <- survfit(Surv(t, event) ~1,
                  data = ipd.b.g)
p.km.b.g <- km_plot(km.b.g)
ggsave("figs/ttd-bsrbr-ipdsurv-eular-good-km.pdf",
       p.km.b.g, height = 5, width = 7)

## corona registry (all patients)
km.c <- survfit(Surv(t, event) ~1, data = ipd.c)
p.km.c <- km_plot(km.c)
ggsave("figs/ttd-corrona-ipdsurv-km.pdf",
       p.km.c, height = 5, width = 7)

# ADJUST BSRBR SURVIVAL CURVES -------------------------------------------------
# cumulative hazard functions
spline.c <- surv_splines(ipd.c, time = seq(0, 120, 1))
spline.b.m <- surv_splines(ipd.b.m, time = seq(0, 120, 1))
spline.b.g <- surv_splines(ipd.b.g, time = seq(0, 120, 1))

# function to adjust survival curves
adjusted_eular_surv <- function(t = seq(1, 120, 1), 
                                n_m = 5492, n_g = 2417, 
                                haz_b_m, haz_b_g, haz_c){
  ## set up
  n <- n_m + n_g
  step <- t[2] - t[1]

  ## estimate ratio of hazards and adjust eular specific hazard functions accordingly
  n.c <- length(haz_c)
  haz.b <- (n_m * haz_b_m[1:n.c] +
               n_g * haz_b_g[1:n.c])/n
  ratio <- haz_c/haz.b
  n.ext <- length(t) - length(haz_c)
  ratio <- c(ratio, rep(ratio[n.c], n.ext))
  haz.b.m.adj <- haz_b_m * ratio
  haz.b.g.adj <- haz_b_g * ratio
  
  ## create eular specific adjusted cumulative hazards
  cumhaz.b.m.adj <- cumsum(haz.b.m.adj * step) 
  cumhaz.b.g.adj <- cumsum(haz.b.g.adj * step) 

  ## return individual patient data
  surv.b.m.adj <- data.frame(tk = c(0, t), sk = c(1, exp(-cumhaz.b.m.adj)))
  surv.b.g.adj <- data.frame(tk = c(0, t), sk = c(1, exp(-cumhaz.b.g.adj)))
  ipd.b.m.adj <- reconstruct_ipd(surv.b.m.adj, risk.b.m, 1)
  ipd.b.g.adj <- reconstruct_ipd(surv.b.g.adj, risk.b.g, 1) 
  ipd.b.m.adj$eular <- "Moderate"
  ipd.b.g.adj$eular <- "Good"
  return(data.table(rbind(ipd.b.m.adj, ipd.b.g.adj)))
}

# adjusted data
ipd.b.adj <- adjusted_eular_surv(haz_b_m = spline.b.m$haz,
                               haz_b_g = spline.b.g$haz,
                               haz_c = spline.c$haz)

# kaplan meier curves (compare adjusted to original)
## moderate responders
km.b.m.adj <- survfit(Surv(t, event) ~1,
                data = ipd.b.adj[eular == "Moderate"])
dat <- data.table(time = c(km.b.m$time/12, km.b.m.adj$time/12),
                  surv = c(km.b.m$surv, km.b.m.adj$surv),
                  lab = c(rep("BSRBR", length(km.b.m$time)),
                          rep("BSRBR (Corrona adjusted)", length(km.b.m.adj$time))))
p.km.b.m.adj <- ggplot(dat, aes(x = time, y = surv, col = lab)) + 
  geom_line() + xlab("Year") + ylab("Survival") + scale_color_discrete("") +
  theme(legend.position = "bottom")
ggsave("figs/ttd-bsrbr-ipdsurv-eular-moderate-km-adjusted.pdf",
       p.km.b.m.adj, height = 5, width = 7)

## good responders
km.b.g.adj <- survfit(Surv(t, event) ~1,
                      data = ipd.b.adj[eular == "Good"])
dat <- data.table(time = c(km.b.g$time/12, km.b.g.adj$time/12),
                  surv = c(km.b.g$surv, km.b.g.adj$surv),
                  lab = c(rep("BSRBR", length(km.b.g$time)),
                          rep("BSRBR (Corrona adjusted)", length(km.b.g.adj$time))))
p.km.b.g.adj <- ggplot(dat, aes(x = time/12, y = surv, col = lab)) + 
  geom_line() + xlab("Year") + ylab("Survival") + scale_color_discrete("") +
  theme(legend.position = "bottom")
ggsave("figs/ttd-bsrbr-ipdsurv-eular-good-km-adjusted.pdf",
       p.km.b.g.adj, height = 5, width = 7)

# ADJUST CORRONA SURVIVAL CURVES BY DISEASE ACTIVITY ---------------------------
corrona_by_da_surv <- function(t = seq(1, 120, 1), km, or){
  km.surv.low <- 1 - or2newprob(1 - km$surv,  or[1])
  km.surv.moderate <- km$surv
  km.surv.high <- 1 - or2newprob(1 - km$surv,  or[3])
  km.df <- data.frame(tk = rep(km$time, 3), 
                   sk = c(km.surv.low, km.surv.moderate, km.surv.high),
                   da = rep(c("Low/remission", "Moderate", "High"), 
                                          each = length(km$surv)))
  risk <- data.frame(trisk = 0, nrisk = 3584)
  ipd.low <- reconstruct_ipd(km.df[km.df$da == "Low/remission", c("tk", "sk")], risk, "Low/remission") 
  ipd.moderate <- reconstruct_ipd(km.df[km.df$da == "Moderate", c("tk", "sk")], risk, "Moderate") 
  ipd.high <- reconstruct_ipd(km.df[km.df$da == "High", c("tk", "sk")], risk, "High") 
  ipd <- rbind(ipd.low, ipd.moderate, ipd.high)
  return(list(ipd = data.table(ipd), km = km.df))
}

# PARAMETRIC FITS --------------------------------------------------------------
mods <- c("exponential", "weibull", "gompertz", "gamma", 
          "llogis", "lnorm", "gengamma")
modlabs <- c("Exponential", "Weibull", "Gompertz", "Gamma",
             "Log-logistic", "Lognormal", "Generalized gamma")

# CORRONA
## overall
fits.c <- parametric_fit(mods, ipd.c)
ic.c <- ic_mat(fits.c, modlabs)
fits.c.diag <- parametric_fits_diag(km.c, fits.c, modlabs)
survfit.data.c <- survfit_data(km.c, fits.c, modlabs, long = FALSE)
surv.gg.c <-  summary(fits.c[["gengamma"]], type = "surv", t = seq(1, 36))[[1]]
p.survfits.c <- ggplot(survfit.data.c[mod %in% c("Kaplan-Meier", "Exponential",
                                                 "Generalized gamma")],
                        aes(x = time/12, y = surv, col = mod)) + geom_line() +
  xlab("Years") + ylab("Survival") + scale_colour_discrete("") +
  theme(legend.position="bottom")
p.gengamma.c <- ggplot(survfit.data.c[mod %in% c("Kaplan-Meier","Generalized gamma")],
                       aes(x = time/12, y = surv, col = mod)) + geom_line() +
  xlab("Years") + ylab("Survival") + scale_colour_discrete("") +
  theme(legend.position="bottom")
write.csv(surv.gg.c, file = "tables/ttd-corrona-gengamma-surv.csv",
          row.names = FALSE)
write.csv(ic.c, "tables/ttd-corrona-ipdsurv-ic.csv")
ggsave("figs/ttd-corrona-ipdsurv-gengamma.pdf",
       p.gengamma.c, height = 5, width = 7)
ggsave("figs/ttd-corrona-ipdsurv-parametric.pdf",
       fits.c.diag$p.survfits, height = 7, width = 10)
ggsave("figs/ttd-corrona-ipdcumhaz-parametric.pdf",
       fits.c.diag$p.cumhazfits, height = 7, width = 10)
ggsave("figs/ttd-corrona-ipdsurv-gengamma-exponential.pdf", p.survfits.c,
       height = 5, width = 7)

## by disease activity
ipd.c.da <- corrona_by_da_surv(km = km.c, or = da.or)
fits.c.da <- parametric_fit(mods, ipd.c.da$ipd, xvars = "arm")
ic.c.da <- ic_mat(fits.c.da, modlabs)

### remission/low vs moderate vs high
tseq = seq(0, 120, 1)
newdat <- data.table(arm = "Low/remission")
survit.c.da.low <- summary(fits.c.da[["gengamma"]], t = tseq, newdata = newdat)[[1]]
newdat <- data.table(arm = "Moderate")
survit.c.da.mod <- summary(fits.c.da[["gengamma"]], t = tseq, newdata = newdat)[[1]]
survfit.c.da <-  rbind(data.table(survit.c.da.low, da = "Low/remission"),
                       data.table(survit.c.da.mod, da = "Moderate"))
p.survfit.c.da <- ggplot(survfit.c.da, aes(x = time/12, y = est, col = da)) + 
  geom_line() +  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = da, linetype = NA), alpha = 0.5) +
  xlab("Year since initiating treatment") + ylab("Proportion surviving") +
  theme(legend.position = "bottom") +
  scale_fill_discrete(name = "Disease activity") +
  scale_colour_discrete(guide = FALSE)
ggsave("figs/ttd-corrona-ipdsurv-gengamma-by-da.pdf", p.survfit.c.da,
       height = 5, width = 7)

# BSRBR
## unadjusted
### moderate 
fits.b.m <- parametric_fit(mods, ipd.b.m)
ic.b.m <- ic_mat(fits.b.m, modlabs)
fits.b.m.diag <- parametric_fits_diag(km.b.m, fits.b.m, modlabs)
write.csv(ic.b.m, "tables/ttd-bsrbr-ipdsurv-eular-moderate-ic.csv")
ggsave("figs/ttd-bsrbr-ipdsurv-eular-moderate-gengamma.pdf",
       fits.b.m.diag$p.survfits.gg, height = 5, width = 7)
ggsave("figs/ttd-bsrbr-ipdsurv-eular-moderate-parametric.pdf",
       fits.b.m.diag$p.survfits, height = 7, width = 10)
ggsave("figs/ttd-bsrbr-ipdcumhaz-eular-moderate-parametric.pdf",
       fits.b.m.diag$p.cumhazfits, height = 7, width = 10)

### good 
fits.b.g <- parametric_fit(mods, ipd.b.g)
ic.b.g <- ic_mat(fits.b.g, modlabs)
fits.b.g.diag <- parametric_fits_diag(km.b.g, fits.b.g, modlabs)
write.csv(ic.b.m, "tables/ttd-bsrbr-ipdsurv-eular-moderate-ic.csv")
ggsave("figs/ttd-bsrbr-ipdsurv-eular-good-gengamma.pdf",
       fits.b.m.diag$p.survfits.gg, height = 5, width = 7)
ggsave("figs/ttd-bsrbr-ipdsurv-eular-good-parametric.pdf",
       fits.b.m.diag$p.survfits, height = 7, width = 10)
ggsave("figs/ttd-bsrbr-ipdcumhaz-eular-good-parametric.pdf",
       fits.b.m.diag$p.cumhazfits, height = 7, width = 10)

### overall
fits.b <- parametric_fit(mods, ipd.b, xvars = "eular")
ic.b <- ic_mat(fits.b, modlabs)

### moderate vs. good 
ic.b <- cbind(ic.b.m, ic.b.g)
colnames(ic.b) <- c("AIC (moderate)", "BIC (moderate)", "AIC (good)", "BIC (good)")
write.csv(ic.b, "tables/ttd-bsrbr-ipdsurv-eular-ic.csv")

survfit.eular <- survfit_by_eular(fits.b.m, fits.b.g, t = seq(0, 120, 1),
                                        dist = "gengamma")
ggsave("figs/ttd-bsrbr-ipdsurv-gengamma-by-eular.pdf",
       survfit.eular$p, height = 7, width = 10)

## adjusted 
### moderate 
fits.b.m.adj <- parametric_fit(mods, ipd.b.adj[eular == "Moderate"])
ic.b.m.adj <- ic_mat(fits.b.m.adj, modlabs)
fits.b.m.adj.diag <- parametric_fits_diag(km.b.m.adj, fits.b.m.adj, modlabs)
write.csv(ic.b.m.adj, "tables/ttd-bsrbr-ipdsurv-eular-moderate-adjusted-ic.csv")
ggsave("figs/ttd-bsrbr-ipdsurv-eular-moderate-gengamma-adjusted.pdf",
       fits.b.m.adj.diag$p.survfits.gg, height = 5, width = 7)
ggsave("figs/ttd-bsrbr-ipdsurv-eular-moderate-adjusted-parametric.pdf",
       fits.b.m.adj.diag$p.survfits, height = 7, width = 10)
ggsave("figs/ttd-bsrbr-ipdcumhaz-eular-moderate-adjusted-parametric.pdf",
       fits.b.m.adj.diag$p.cumhazfits, height = 7, width = 10)

### good
fits.b.g.adj <- parametric_fit(mods, ipd.b.adj[eular == "Good"])
ic.b.g.adj <- ic_mat(fits.b.g.adj, modlabs)
fits.b.g.adj.diag <- parametric_fits_diag(km.b.g.adj, fits.b.g.adj, modlabs)
write.csv(ic.b.g.adj, "tables/ttd-bsrbr-ipdsurv-eular-good-adjusted-ic.csv")
ggsave("figs/ttd-bsrbr-ipdsurv-eular-good-gengamma-adjusted.pdf",
       fits.b.g.adj.diag$p.survfits.gg, height = 5, width = 7)
ggsave("figs/ttd-bsrbr-ipdsurv-eular-good-adjusted-parametric.pdf",
       fits.b.g.adj.diag$p.survfits, height = 7, width = 10)
ggsave("figs/ttd-bsrbr-ipdcumhaz-eular-good-adjusted-parametric.pdf",
       fits.b.g.adj.diag$p.cumhazfits, height = 7, width = 10)

### moderate vs. good
ic.b.adj <- cbind(ic.b.m.adj, ic.b.g.adj)
colnames(ic.b.adj) <- c("AIC (moderate)", "BIC (moderate)", "AIC (good)", "BIC (good)")
write.csv(ic.b.adj, "tables/ttd-bsrbr-ipdsurv-eular-adjusted-ic.csv")

survfit.eular.adj <- survfit_by_eular(fits.b.m.adj, fits.b.g.adj, t = seq(0, 120, 1),
                                        dist = "gengamma")
ggsave("figs/ttd-bsrbr-ipdsurv-adjusted-gengamma-by-eular.pdf",
       survfit.eular.adj$p, height = 7, width = 10)

## adjusted and unadjusted
survfit.eular.comp <- rbind(data.table(survfit.eular$data, lab = "BSRBR"),
                            data.table(survfit.eular.adj$data, lab = "BSRBR (Corrona adjusted)"))
p.survfit.eular.comp <- ggplot(survfit.eular.comp, aes(x = time/12, y = est, col = eular)) + 
  geom_line() +  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = eular, linetype = NA), alpha = 0.5) +
  facet_wrap(~lab) +   xlab("Year since initiating treatment") + ylab("Proportion surviving") +
  theme(legend.position = "bottom") +
  scale_fill_discrete(name = "Eular response") +
  scale_colour_discrete(guide = FALSE)
ggsave("figs/ttd-bsrbr-ipdsurv-comp-gengamma-by-eular.pdf",
       p.survfit.eular.comp, height = 7, width = 10)

# SAVE PARAMETERS FOR MODEL ----------------------------------------------------
# moderate responders
ttd.eular.mod <- ttd.eular.good <- ttd.eular.mod.adj <- ttd.eular.good.adj <- list()
ttd.da <- ttd.all <- list()
for (i in 1:length(mods)){
  ttd.eular.mod[[mods[i]]] <- flexsurvreg_pars(fits.b.m[[mods[i]]])
  ttd.eular.good[[mods[i]]] <- flexsurvreg_pars(fits.b.g[[mods[i]]])
  ttd.eular.mod.adj[[mods[i]]] <- flexsurvreg_pars(fits.b.m.adj[[mods[i]]])
  ttd.eular.good.adj[[mods[i]]] <- flexsurvreg_pars(fits.b.g.adj[[mods[i]]])
  ttd.all[[mods[i]]] <- flexsurvreg_pars(fits.c[[mods[i]]])
  ttd.da[[mods[i]]] <- flexsurvreg_pars(fits.c.da[[mods[i]]])
}
names(ttd.all$exponential$est) <- "rate"
names(ttd.eular.mod$exponential$est) <- "rate"
names(ttd.eular.good$exponential$est) <- "rate"
ttd.eular.uk <- list(moderate = ttd.eular.mod, good = ttd.eular.good)
ttd.eular <- list(moderate = ttd.eular.mod.adj, good = ttd.eular.good.adj) # us version

# SAVE PARAMETERS --------------------------------------------------------------
save(ttd.all, ttd.eular, ttd.da,
     file = "../data/ttd.rda", compress = "bzip2")
