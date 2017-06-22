rm(list = ls())
eular2haq <- data.frame(mean = c(0, -.317, -.672),
                               se = c(0, .048, .112))
save(eular2haq, 
     file = "../data/eular2haq.rda", compress = "bzip2")