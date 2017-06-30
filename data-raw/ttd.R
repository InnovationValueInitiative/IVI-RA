rm(list = ls())
library("survival")
library("flexsurv")
library("ggplot2")
library("MASS")
library("data.table")
source("func.R")
theme_set(theme_bw())

# LOAD DIGITIZED SURVIVAL CURVES AND NUMBERS AT RISK ---------------------------
surv.b.g <- read.csv("ttd-bsrbr-eular-good-ggamma.csv")
risk.b.g <- read.csv("ttd-bsrbr-eular-good-risk.csv")
surv.b.m <- read.csv("ttd-bsrbr-eular-moderate-ggamma.csv")
risk.b.m <- read.csv("ttd-bsrbr-eular-moderate-risk.csv")
surv.c <- read.csv("ttd-corrona-allpat.csv")
risk.c <- read.csv("ttd-corrona-allpat-risk.csv")

# CONVERT CURVES TO INDIVIDUAL PATIENT LEVEL PATIENT DATA ----------------------
## algorithm
ipd.b.g <- reconstruct_ipd(surv.b.g, risk.b.g, 1)
ipd.b.g$t <- (ipd.b.g$t/30.42) * 1000 # convert time from 1000s of days to months
ipd.b.m <- reconstruct_ipd(surv.b.m, risk.b.m, 1)
ipd.b.m$t <- (ipd.b.m$t/30.42) * 1000
ipd.c <- reconstruct_ipd(surv.c, risk.c, 1, tot_events = 3584)

## check kaplan meier of curves
# moderate responders
km.b.m <- survfit(Surv(t, event) ~1,
                       data = ipd.b.m)
p.km.b.m <- km_plot(km.b.m)
ggsave("figs/ttd-bsrbr-ipdsurv-eular-moderate-km.pdf",
       p.km.b.m, height = 5, width = 7)

# good responders
km.b.g <- survfit(Surv(t, event) ~1,
                  data = ipd.b.g)
p.km.b.g <- km_plot(km.b.g)
ggsave("figs/ttd-bsrbr-ipdsurv-eular-good-km.pdf",
       p.km.b.g, height = 5, width = 7)

# corona registry (all patients)
km.c <- survfit(Surv(t, event) ~1,
                     data = ipd.c)
p.km.c <- km_plot(km.c)
ggsave("figs/ttd-bsrbr-ipdsurv-corrona-km.pdf",
       p.km.c, height = 5, width = 7)

# ADJUSTED SURVIVAL CURVES -----------------------------------------------------
## cumulative hazard functions
spline.c <- surv_splines(ipd.c, time = seq(0, 120, 1))
spline.b.m <- surv_splines(ipd.b.m, time = seq(0, 120, 1))
spline.b.g <- surv_splines(ipd.b.g, time = seq(0, 120, 1))

## function to adjust survival curves
adjusted_eular_surv <- function(t = seq(1, 120, 1), 
                                n_m = 5492, n_g = 2417, 
                                haz_b_m, haz_b_g, haz_c){
  # set up
  n <- n_m + n_g
  step <- t[2] - t[1]

  # estimate ratio of hazards and adjust eular specific hazard functions accordingly
  n.c <- length(haz_c)
  haz.b <- (n_m * haz_b_m[1:n.c] +
               n_g * haz_b_g[1:n.c])/n
  ratio <- haz_c/haz.b
  n.ext <- length(t) - length(haz_c)
  ratio <- c(ratio, rep(ratio[n.c], n.ext))
  haz.b.m.adj <- haz_b_m * ratio
  haz.b.g.adj <- haz_b_g * ratio
  
  # create eular specific adjusted cumulative hazards
  cumhaz.b.m.adj <- cumsum(haz.b.m.adj * step) 
  cumhaz.b.g.adj <- cumsum(haz.b.g.adj * step) 

  # return individual patient data
  surv.b.m.adj <- data.frame(tk = c(0, t), sk = c(1, exp(-cumhaz.b.m.adj)))
  surv.b.g.adj <- data.frame(tk = c(0, t), sk = c(1, exp(-cumhaz.b.g.adj)))
  ipd.b.m.adj <- reconstruct_ipd(surv.b.m.adj, risk.b.m, 1)
  ipd.b.g.adj <- reconstruct_ipd(surv.b.g.adj, risk.b.g, 1) 
  return(list(moderate = ipd.b.m.adj, good = ipd.b.g.adj))
}

## adjusted data
ipd.b.adj <- adjusted_eular_surv(haz_b_m = spline.b.m$haz,
                               haz_b_g = spline.b.g$haz,
                               haz_c = spline.c$haz)

## kaplan meier curves (compare adjusted to original)
# moderate responders
km.b.m.adj <- survfit(Surv(t, event) ~1,
                data = ipd.b.adj$moderate)
dat <- data.table(time = c(km.b.m$time/12, km.b.m.adj$time/12),
                  surv = c(km.b.m$surv, km.b.m.adj$surv),
                  lab = c(rep("BSRBR", length(km.b.m$time)),
                          rep("BSRBR (Corrona adjusted)", length(km.b.m.adj$time))))
p.km.b.m.adj <- ggplot(dat, aes(x = time, y = surv, col = lab)) + 
  geom_line() + xlab("Year") + ylab("Survival") + scale_color_discrete("") +
  theme(legend.position = "bottom")
ggsave("figs/ttd-bsrbr-ipdsurv-eular-moderate-km-adjusted.pdf",
       p.km.b.m.adj, height = 5, width = 7)

# good responders
km.b.g.adj <- survfit(Surv(t, event) ~1,
                      data = ipd.b.adj$good)
dat <- data.table(time = c(km.b.g$time/12, km.b.g.adj$time/12),
                  surv = c(km.b.g$surv, km.b.g.adj$surv),
                  lab = c(rep("BSRBR", length(km.b.g$time)),
                          rep("BSRBR (Corrona adjusted)", length(km.b.g.adj$time))))
p.km.b.g.adj <- ggplot(dat, aes(x = time/12, y = surv, col = lab)) + 
  geom_line() + xlab("Year") + ylab("Survival") + scale_color_discrete("") +
  theme(legend.position = "bottom")
ggsave("figs/ttd-bsrbr-ipdsurv-eular-good-km-adjusted.pdf",
       p.km.b.g.adj, height = 5, width = 7)

# PARAMETRIC FITS --------------------------------------------------------------
mods <- c("exponential", "weibull", "gompertz", "gamma", 
          "llogis", "lnorm", "gengamma")
modlabs <- c("Exponential", "Weibull", "Gompertz", "Gamma",
             "Log-logistic", "Lognormal", "Generalized gamma")

## unadjusted
# moderate 
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

# good 
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

# moderate vs. good 
ic.b <- cbind(ic.b.m, ic.b.g)
colnames(ic.b) <- c("AIC (moderate)", "BIC (moderate)", "AIC (good)", "BIC (good)")
write.csv(ic.b, "tables/ttd-bsrbr-ipdsurv-eular-ic.csv")

survfit.eular <- survfit_by_eular(fits.b.m, fits.b.g, t = seq(0, 120, 1),
                                        dist = "gengamma")
ggsave("figs/ttd-bsrbr-ipdsurv-gengamma-by-eular.pdf",
       survfit.eular$p, height = 7, width = 10)

## adjusted 
# moderate 
fits.b.m.adj <- parametric_fit(mods, ipd.b.adj$moderate)
ic.b.m.adj <- ic_mat(fits.b.m.adj, modlabs)
fits.b.m.adj.diag <- parametric_fits_diag(km.b.m.adj, fits.b.m.adj, modlabs)
write.csv(ic.b.m.adj, "tables/ttd-bsrbr-ipdsurv-eular-moderate-adjusted-ic.csv")
ggsave("figs/ttd-bsrbr-ipdsurv-eular-moderate-gengamma-adjusted.pdf",
       fits.b.m.adj.diag$p.survfits.gg, height = 5, width = 7)
ggsave("figs/ttd-bsrbr-ipdsurv-eular-moderate-adjusted-parametric.pdf",
       fits.b.m.adj.diag$p.survfits, height = 7, width = 10)
ggsave("figs/ttd-bsrbr-ipdcumhaz-eular-moderate-adjusted-parametric.pdf",
       fits.b.m.adj.diag$p.cumhazfits, height = 7, width = 10)

# good
fits.b.g.adj <- parametric_fit(mods, ipd.b.adj$good)
ic.b.g.adj <- ic_mat(fits.b.g.adj, modlabs)
fits.b.g.adj.diag <- parametric_fits_diag(km.b.g.adj, fits.b.g.adj, modlabs)
write.csv(ic.b.g.adj, "tables/ttd-bsrbr-ipdsurv-eular-good-adjusted-ic.csv")
ggsave("figs/ttd-bsrbr-ipdsurv-eular-good-gengamma-adjusted.pdf",
       fits.b.g.adj.diag$p.survfits.gg, height = 5, width = 7)
ggsave("figs/ttd-bsrbr-ipdsurv-eular-good-adjusted-parametric.pdf",
       fits.b.g.adj.diag$p.survfits, height = 7, width = 10)
ggsave("figs/ttd-bsrbr-ipdcumhaz-eular-good-adjusted-parametric.pdf",
       fits.b.g.adj.diag$p.cumhazfits, height = 7, width = 10)

# moderate vs. good
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
for (i in 1:length(mods)){
  ttd.eular.mod[[mods[i]]] <- flexsurvreg_pars(fits.b.m[[mods[i]]])
  ttd.eular.good[[mods[i]]] <- flexsurvreg_pars(fits.b.g[[mods[i]]])
  ttd.eular.mod.adj[[mods[i]]] <- flexsurvreg_pars(fits.b.m.adj[[mods[i]]])
  ttd.eular.good.adj[[mods[i]]] <- flexsurvreg_pars(fits.b.g.adj[[mods[i]]])
}

# SAVE PARAMETERS --------------------------------------------------------------
save(ttd.eular.mod, ttd.eular.good, ttd.eular.mod.adj, ttd.eular.good.adj, 
     file = "../data/ttd-pars.rda", compress = "bzip2")
