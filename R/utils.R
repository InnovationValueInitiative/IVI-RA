#' Convert list of matrices to an array
#'
#' Convert a list of matrices to an array.
#' 
#' @param l List of matrices
#' @export
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

#' ACR response from NMA
#'
#' Calculate ACR response probability from ordered probit NMA
#' 
#' @param A Mean probability that ACR response < 20 from cDMARDS 
#' @param z2 Mean for ACR 50 cutpoint
#' @param z3 Mean for ACR 70 cutpoint
#' @param delta Vector of regression coefficient means (on probit scale) for therapies
#' @export
acr_prob <- function(A, z2, z3, delta){
  p <- matrix(NA, nrow = length(delta), ncol = 4)
  pl <- rep(NA, 3)
  for (i in 1:length(delta)){
    # probability less than category
    pl[3] <- pnorm(A + z3 + delta[i])
    pl[2] <- pnorm(A + z2 + delta[i])
    pl[1] <- pnorm(A + delta[i])
    
    # probability in category
    p[i, 1] <- pl[1]
    p[i, 2] <- pl[2] - pl[1]
    p[i, 3] <- pl[3] - pl[2]
    p[i, 4] <- 1 - pl[3]
  }
  rownames(p) <- names(delta)
  p[which(row.names(p) == "placebo"), ] <- c(1, 0, 0, 0)
  p[which(row.names(p) == "nbt"), ] <- c(1, 0, 0, 0)
  return(p)
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