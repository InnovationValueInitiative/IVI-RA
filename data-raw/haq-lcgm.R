rm(list = ls())
library("ggplot2")
library("data.table")
theme_set(theme_bw())

# NORTON (2014) LATENT CLASS GROWTH MODEL (LCGM) ------------------------------
# predicted HAQ
haq <- data.table(read.table("haqpred-norton-2014.dat"))
setnames(haq, colnames(haq),
         c("time", "c3_obs", "c3_exp", "c1_obs", "c1_exp", 
           "c4_obs", "c4_exp", "c2_obs", "c2_exp"))
haq[, time := c(0, .5, seq(1, 15))]
haq[, xt := 1 - (1/(time + 1))]
c1.lm <- lm(c1_obs ~ xt + I(xt^2) + I(xt^3), data = haq[time >= 2])
c2.lm <-  lm(c2_obs ~ xt + I(xt^2) + I(xt^3), data = haq[time >= 2, ])
c3.lm <-  lm(c3_obs ~ xt + I(xt^2) + I(xt^3), data = haq[time >= 2, ])
c4.lm <-  lm(c4_obs ~ xt + I(xt^2) + I(xt^3), data = haq[time >= 2, ])
haq[, c1_lmexp := c(rep(NA, 3), predict(c1.lm))]
haq[, c2_lmexp := c(rep(NA, 3), predict(c2.lm))]
haq[, c3_lmexp := c(rep(NA, 3), predict(c3.lm))]
haq[, c4_lmexp := c(rep(NA, 3), predict(c4.lm))]
p.dat <- melt(haq[, .(time, c1_obs, c2_obs, c3_obs, c4_obs, c1_lmexp, c2_lmexp,
                      c3_lmexp, c4_lmexp)], id.vars = "time")
p.dat[, class := ifelse(variable %in% c("c1_obs", "c1_lmexp"), 1, NA)]
p.dat[, class := ifelse(variable %in% c("c2_obs", "c2_lmexp"), 2, class)]
p.dat[, class := ifelse(variable %in% c("c3_obs", "c3_lmexp"), 3, class)]
p.dat[, class := ifelse(variable %in% c("c4_obs", "c4_lmexp"), 4, class)]
p.dat[, type := ifelse(grepl("lmexp", variable) == TRUE, "Expected", "Observed")]
p <- ggplot(p.dat, aes(x = time, y = value, linetype = factor(type), col = factor(class))) + 
  geom_line() +
  geom_point(data = p.dat[ type == "Observed"], aes(x = time, y = value), size = 1) +
  xlab("Follow-up in years") + ylab("HAQ") + scale_linetype_discrete("Value") +
  scale_color_discrete("Class") +
  theme(legend.position = "bottom")
ggsave("figs/haq-lcgm-obsexp.pdf", p, height = 5, width = 7)

# latent class probabilities
coefs <- fread("coef-norton-2014.csv")
beta1 <- coef(c1.lm)
beta1.se <- sqrt(diag(vcov(c1.lm)))
beta2 <- coef(c2.lm)
beta2.se <- sqrt(diag(vcov(c2.lm)))
beta3 <- coef(c3.lm)
beta3.se <- sqrt(diag(vcov(c3.lm)))
beta4 <- coef(c4.lm)
beta4.se <- sqrt(diag(vcov(c4.lm)))
beta.var <- paste0(c("intercept_", "linear_", "quadratic_", "cubic_"), rep(seq(1,4), each = 4))
beta.parameter <- rep(c("beta1", "beta2", "beta3", "beta4"), each = 4)
beta.dt <- data.table(var = beta.var, data = "eras", parameter = beta.parameter, 
                      coef = c(beta1, beta2, beta3, beta4),
                      se = c(beta1.se, beta2.se, beta3.se, beta4.se))
coefs <- rbind(beta.dt, coefs)
index <- list()
parameters <- unique(coefs$parameter)
for (par in parameters){
  index[[par]] <- which(coefs$parameter == par)
}

# save
haq.lcgm.pars <- list(coef = coefs, index = index)
save(haq.lcgm.pars, 
     file = "../data/haq-lcgm-pars.rda", compress = "bzip2")
