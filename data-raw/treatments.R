rm(list = ls())
library(data.table)
treatments <- fread("treatments.csv")
treatments[, route := factor(route)]
treatments[, approval_date := ifelse(sname %in% c("nbt", "placebo"), 
                                     "12/31/88", approval_date)]
treatments[, years_since_approval := as.numeric(difftime(as.Date(end_date, "%m/%d/%y"), 
                                                         as.Date(approval_date, "%m/%d/%y"), 
                                                         unit = "weeks"))/52.25]
treatments[, c("end_date") := NULL]
save(treatments,
     file = "../data/treatments.rda", compress = "bzip2")

