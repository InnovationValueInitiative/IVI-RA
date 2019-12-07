rm(list = ls())
library("data.table")
source("func.R")
mort.hr.haqdif <- fread("mort-hr-haqdif.csv")
mort.hr.haqdif[, loghr := log(hr)]
mort.hr.haqdif[, loghr_lower := log(hr_lower)]
mort.hr.haqdif[, loghr_upper := log(hr_upper)]
mort.hr.haqdif[, loghr_se := se_normal(loghr_lower, loghr_upper)]
save(mort.hr.haqdif, file = "../data/mort-hr-haqdif.rda", compress = "bzip2")