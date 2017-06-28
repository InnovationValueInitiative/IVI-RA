rm(list = ls())
library("data.table")
acr2haq <- data.table(acr_cat = c("ACR < 20", "ACR 20-50", "ACR 50-70", "ACR 70+" ),
                      mean = c(.11, .44, .76, 1.07),
                      se = c(.06765, .05657, .09059, .07489))
save(acr2haq, 
     file = "../data/acr2haq.rda", compress = "bzip2")
