# clean raw nma output
rm(list = ls())
library("data.table")


# ACR response
acr <- fread("nma-raw/nma-acr-naive-re-coda1.csv")
acr <- acr[, grep('^(A|d\\[|z)', colnames(acr)), with = FALSE]
write.csv(acr, "nma-acr-naive-re-coda.csv", row.names = FALSE)

# HAQ change
haq <- fread("nma-raw/nma-haq-naive-re-coda1.csv")
haq <- haq[, grep('^(A|d\\[|z)', colnames(haq)), with = FALSE]
write.csv(haq, "nma-haq-naive-re-coda.csv", row.names = FALSE)

# DAS28 change
das28 <- fread("nma-raw/nma-das28-naive-re-coda1.csv")
das28 <- das28[, grep('^(A|d\\[|z)', colnames(das28)), with = FALSE]
write.csv(das28, "nma-das28-naive-re-coda.csv", row.names = FALSE)
