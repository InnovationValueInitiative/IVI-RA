age <- 55.2
female <- .7
das28 <- 5.46
disease_duration <- 8.2/12
rf <- .73
acr <- .706
ses <- .186

e1 <- 1
e2 <- exp(-3.496 + .025 * age + .841 * female + 0.304 * das28 + .032 * disease_duration + 
            .214 * rf + .278 * acr + .993 * ses)
e3 <- exp(-6.686 + 0.037 * age + 1.694 * female + 0.573 * das28 + 0.046 * disease_duration +
            0.315 * rf + 0.413 * acr + 1.119 * ses)
e4 <- exp(-12.055 + 0.082 * age + 1.976 * female + 0.800 * das28 + 0.042 * disease_duration + 
            0.298* rf + 0.939 * acr + 1.429 * ses)
p1 <- e1/(e1 + e2 + e3 + e4)
p2 <- e2/(e1 + e2 + e3 + e4)
p3 <- e3/(e1 + e2 + e3 + e4)
p4 <- e4/(e1 + e2 + e3 + e4)

print(c(p1, p2, p3, p4))