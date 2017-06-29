rm(list = ls())
library("data.table")
eular2haq <- data.table(eular = c("None", "Moderate", "Good"),
                        mean = c(0, -.317, -.672),
                        se = c(0, .048, .112))
save(eular2haq, 
     file = "../data/eular2haq.rda", compress = "bzip2")