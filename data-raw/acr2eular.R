rm(list = ls())
acr2eular <- matrix(c(755, 4, 2, 0, 136, 27, 2, 2, 57, 26, 10, 2), 
                    nrow = 4, ncol = 3, byrow = FALSE)
rownames(acr2eular) <- c("acr_sub", "acr20", "acr50", "acr70")
colnames(acr2eular) <- c("eular_none", "eular_moderate", "eular_good")
save(acr2eular, 
     file = "../data/acr2eular.rda", compress = "bzip2")