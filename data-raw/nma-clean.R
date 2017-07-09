# clean raw nma output
rm(list = ls())
library("data.table")
acr <- fread("nma-raw/bio naive coda 1 re.csv")
acr <- acr[, grep('^(A|d\\[|z)', colnames(acr)), with = FALSE]
write.csv(acr, "nma-acr-naive-re-coda.csv", row.names = FALSE)