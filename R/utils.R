#' Convert list of matrices to an array
#'
#' Convert a list of matrices to an array.
#' 
#' @param l List of matrices
#' @keywords internal
list2array <- function(l){
  if (class(l[[1]]) == "numeric"){
      dim1 <- 1; dim2 <- length(l[[1]]); dim3 <- length(l)
  } else {
      dim1 <- nrow(l[[1]]); dim2 <- ncol(l[[1]]); dim3 <- length(l) 
  }
  x <-  array(unlist(l),
                dim = c(dim1, dim2, dim3))
 
  return(x)
}

#' NMA parameters to change in continuous outcome
#'
#' Calculate change in continuous outcome from linear model used for NMA.
#' 
#' @param A A vector from the posterior distribution of the probability that ACR response < 20 
#' from cDMARDS.
#' @param delta A matrix from the posterior distribution of regression coefficients 
#' for each therapy in the NMA.
#' 
#' @return Matrix of change in continuous outcome variable for each therapy. 
#' 
#' @export
nma_lm2prob <- function(A, delta){
  nther <- ncol(delta)
  nsims <- length(A)
  dy <- matrix(NA, nrow = nsims, ncol = nther)
  for (i in 1:nther){
    dy[, i] <- A +  delta[, i]
    #dy[, i] <- ifelse(dy[, i] < 0, rr * dy[, i], 1/rr * dy[, i])
  }
  return(dy)
} 

#' NMA parameters to ACR response probabilities
#'
#' Calculate ACR response probabilities from NMA ordered probit parameters. 
#' 
#' @param A A vector from the posterior distribution of the probability that ACR response < 20 
#' from cDMARDS.
#' @param z2 A vector from the posterior distribution for the ACR 50 cutpoint.
#' @param z3 A vector from the posterior distribution for the ACR 70 cutpoint.
#' @param delta A matrix from the posterior distribution of regression coefficients 
#' (on the probit scale) for each therapy in the NMA.
#' @param rr Sampled value of elative risk.
#' 
#' @return List containing an array \code{non.overlap} (the probability of ACR < 20, ACR 20-50,
#' ACR 50-70, and ACR 70+) and \code{overlap} (the probability of ACR20/50/70).
#' 
#' @export
nma_acrprob <- function(A, z2, z3, delta, rr = 1){
  nther <- ncol(delta)
  nsims <- length(A)
  p <- array(NA, dim = c(nsims, 4, nther))
  pl <- matrix(NA, nrow = nsims, ncol = 3)
  po <- array(NA, dim = c(nsims, 4, nther))
  if (is.numeric(pl)) pl <- t(as.matrix(pl))
  for (i in 1:nther){
    # probability less than category
    pl[, 3] <- pnorm(A + z3 + delta[, i]) # less than ACR 70
    pl[, 2] <- pnorm(A + z2 + delta[, i]) # less than ACR 50
    pl[, 1] <- pnorm(A + delta[, i]) # less than ACR 20
    
    # probability in overlapping categories
    po[, 4, i] <- rr * (1 - pl[, 3]) # greater than ACR 70
    po[, 3, i] <- rr * (1 - pl[, 2])  # greater than ACR 50
    po[, 2, i] <- rr * (1 - pl[, 1]) # greater than ACR 20
    po[, 1, i] <- 1 - po[, 2, i] # less than ACR 20
    
    # probability in mutually exclusive categories
    p[, 1, i] <- po[, 1, i] # less than ACR 20
    p[, 2, i] <- po[, 2, i] - po[, 3, i] # ACR 20-50
    p[, 3, i] <- po[, 3, i] - po[, 4, i]  # ACR 50-70
    p[, 4, i] <- po[, 4, i] # ACR 70
  }
  return(list(non.overlap = p, overlap = po))
} 

#' Simulate survival time for RA patients
#'
#' Simulate survival time given age and a vector of HAQ scores
#' 
#' @param n Number of observations
#' @param age Age of patients. Length must be 1.
#' @param haq0 Baseline HAQ. Length must be 1.
#' @param male Vector indicating patient gender (1 = male, 0 = female)
#' @param haq Vector of HAQ scores (e.g. HAQ trajectory)
#' @param ltfemale Lifetable for women. Must contain column 'age' for single-year of age and 'qx' for
#' the probability of death at a given age. 
#' @param ltmale Identical to \code{ltfemale} but for men.
#' @param x_mort Design matrix for mortality adjustment with odds ratios. Number of rows equal to n. 
#' @param logor Vector of log odds ratio mortality adjusters. Length equal to number of columns in x.
#' @param cycle_length Length of model cycles in months. Default is 1 month cycles.
#' @param loghr Log hazard ratio of impact of change in HAQ from baseline on mortality rate. A vector with
#' each element denoting (in order) hazard ratio for months 0-6, >6 - 12, >12 - 24, >24 -36, >36.
#' @export
sim_survtime <- function(n = 1000, age = 55, haq0 = 1, male = 0, 
                      haq = rep(haq0, 12/cycle_length * (100 - age + 1)), 
                      ltfemale = lifetable.female, ltmale = lifetable.male, 
                      x_mort = haq0, logor = mort.or$logor, 
                      cycle_length = 6, loghr = mort.hr.haqdif$loghr){
  lt.data <- lt_data(ltfemale, ltmale)
  le <- sample_survC(n, age, male, lt.data$male, lt.data$female,
                   x_mort, logor, haq0, haq,
                   cycle_length, loghr)
  return(le)
}