context("Mortality functions")

# Test C++ function newprobC --------------------------------------------------
inv_logit <- function(x){
  return(exp(x)/(1 + exp(x)))
}

test_that("newprob", {
  x <- c(1, 1, 0)
  logor <- c(1.1, 1.3, .85)
  logit.baseline <- .5
  predR <- as.numeric(inv_logit(x %*%  logor + logit.baseline))
  predC <- IVI026:::newprobC(x, logor, logit.baseline)
  expect_equal(predR, predC)
})

# Test C++ function mortprobC -------------------------------------------------
# parameters
bhaq <- 1.5
chaq <- 2.0
x.mort <- as.matrix(bhaq)
mort.logor <- .79
qx <- .01
cycle.len <- 6
age <- 55
loghr <- mort.hr.haqdif$loghr
lt <- lt_data(ltfemale = lifetable.female, ltmale = lifetable.male)

# simple R function
mortprobR <- function(age, male, lifetable_male, lifetable_female,                      
                      x, logor, haq0, haq, cycle_length, month, loghr){
  if (age >= 100){
      qx <- 1
  } else {
      if (male == 1){
          agerow <- which(lifetable_male[, "age"] == 55)
          logit.qx <- lifetable_male[agerow, "logit_qx"]
      } else{
          agerow <- which(lifetable_female[, "age"] == 55)
          logit.qx <- lifetable_female[agerow, "logit_qx"]
      }
      qx <- IVI026:::newprobC(x, logor, logit.qx)
      haq.change <- (haq - haq0)/.25
      rate <- -log(1 - qx) * exp(loghr * haq.change) 
      qx <- 1 - exp(-rate * (cycle_length/12))
  }
  return(qx)
}

# test
mort.pars <- list(age, male = 0, lifetable_male = lt$male,
                   lifetable_female = lt$female,
                   x = x.mort, logor = mort.logor, haq0 = bhaq, haq = chaq,
                   cycle_length = cycle.len, month = 4, loghr)
mort.parsR <- mort.pars
mort.parsR[[length(mort.parsR)]] <- mort.parsR[[length(mort.parsR)]][1]
test_that("mortprobC", {
  qxC <- do.call(getFromNamespace("mortprobC", "IVI026"), mort.pars)
  qxR <- do.call("mortprobR", mort.pars)
  expect_equal(qxC, qxR[1])
})

# Test C++ function death_sampleC ---------------------------------------------
death_sampleR <- function(n, ...){
  qx <- mortprobR(...)
  return(rbinom(n, 1, qx))
}

test_that("death_sampleC", {
  # C++ function
  set.seed(100)
  sampC <- replicate(1000, do.call(getFromNamespace("sample_deathC", "IVI026"), mort.pars))
  
  # R function
  set.seed(100)
  sampR <- do.call("death_sampleR", c(1000, mort.parsR))
  expect_equal(sampC, sampR)
})