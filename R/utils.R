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
#' @param d A matrix from the posterior distribution of the "d's" 
#' for each treatment in the NMA.
#' 
#' @return Matrix of change in continuous outcome variable for each therapy. 
#' 
#' @export
nma_lm_dy <- function(A, d){
  nther <- ncol(d)
  nsims <- length(A)
  dy <- matrix(NA, nrow = nsims, ncol = nther)
  for (i in 1:nther){
    dy[, i] <- A +  d[, i]
    #dy[, i] <- ifelse(dy[, i] < 0, rr * dy[, i], 1/rr * dy[, i])
  }
  return(dy)
} 

#' ACR response probabilities
#'
#' Calculate ACR response probabilities from an ordered probit model of ACR response used for a
#' network meta-analysis. 
#' 
#' @param A A vector from the posterior distribution of the probability that ACR response < 20 
#' from cDMARDS.
#' @param z2 A vector from the posterior distribution for the ACR 50 cutpoint.
#' @param z3 A vector from the posterior distribution for the ACR 70 cutpoint.
#' @param d A matrix from the posterior distribution of the "d's" 
#' for each treatment in the NMA.
#' @param k Sampled value of constant \eqn{k}.
#' 
#' @return List containing an array \code{non.overlap} (the probability of ACR < 20, ACR 20-50,
#' ACR 50-70, and ACR 70+) and \code{overlap} (the probability of ACR20/50/70).
#' 
#' @export
nma_acrprob <- function(A, z2, z3, d, rr = 1){
  nther <- ncol(d)
  nsims <- length(A)
  p <- array(NA, dim = c(nsims, 4, nther))
  pl <- matrix(NA, nrow = nsims, ncol = 3)
  po <- array(NA, dim = c(nsims, 4, nther))
  if (is.numeric(pl)) pl <- t(as.matrix(pl))
  for (i in 1:nther){
    # probability less than category
    pl[, 3] <- pnorm(A + z3 + d[, i]) # less than ACR 70
    pl[, 2] <- pnorm(A + z2 + d[, i]) # less than ACR 50
    pl[, 1] <- pnorm(A + d[, i]) # less than ACR 20
    
    # probability in overlapping categories
    po[, 4, i] <- k * (1 - pl[, 3]) # greater than ACR 70
    po[, 3, i] <- k * (1 - pl[, 2])  # greater than ACR 50
    po[, 2, i] <- k * (1 - pl[, 1]) # greater than ACR 20
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
#' each element denoting (in order) hazard ratio for \eqn{month \le 6}, \eqn{6 > month \le 12}, 
#'   \eqn{12 < month \le 24}, \eqn{24 < month \le 36}, \eqn{month > 36}.
#'   
#' @examples 
#' set.seed(101)
#' t <- seq(0, 45 * 2, .5) # 90 model cycles, or 45 years.
#' haq.seq <-  pmin(1.5 + .03 * t, 3)
#' surv <- sim_survtime(n = 1000, age = 55, haq = haq.seq, male = 0)
#' summary(surv)
#'       
#' @export
sim_survtime <- function(n = 1000, age = 55, haq0 = 1, male = 0, 
                      haq = rep(haq0, 12/cycle_length * (100 - age + 1)), 
                      ltfemale = lifetable.female, ltmale = lifetable.male, 
                      x_mort = haq0, logor = mort.or$logor, 
                      cycle_length = 6, loghr = mort.hr.haqdif$loghr){
  lt.data <- lt_pars(ltfemale, ltmale)
  le <- sample_survC(n, age, male, lt.data$male, lt.data$female,
                   x_mort, logor, haq0, haq,
                   cycle_length, loghr)
  return(le)
}

#' Check scalar
#'
#' Check argument that must be a scalar is correct. 
#' 
#' 
#' @return Error message if incorrect arguments are passed.
#'  
#' @keywords internal
check_scalar <- function(s, pos = FALSE){
  s.name <- deparse(substitute(s))
  if(!is.atomic(s) & length(s) == 1L){
    stop(paste0(s.name, " must be a scalar of length 1."))
  }
  if (pos == TRUE){
    if(s < 0){
      stop(paste0(s.name, " must be greater than or equal to 0."))
    }
  }
}

#' Check vector
#'
#' Check argument that must be a vector is correct. 
#' 
#' 
#' @return Error message if incorrect arguments are passed.
#'  
#' @keywords internal
check_vector <- function(v, len = NULL, pos = FALSE){
  v.name <- deparse(substitute(v))
  if(!is.vector(v)){
    stop(paste0(v.name, " must be a vector."))
  }
  if(!is.null(len)){
    if(length(v) != len){
      stop(paste0(v.name, "  must be length ", len, "."))
    }
  }
  if (pos == TRUE){
    if(any(v < 0)){
      stop(paste0("Each element of ", v.name, " must be greater than or equal to 0."))
    }
  }
}

#' Check matrix
#'
#' Check argument that must be a matrix is correct. 
#' 
#' 
#' @return Error message if incorrect arguments are passed.
#'  
#' @keywords internal
check_matrix <- function(m, nrow, ncol, pos = FALSE){
  m.name <- deparse(substitute(m))
  if(!is.matrix(m)){
    stop(paste0(m.name, " must be a matrix"))
  }
  if(nrow(m) != nrow){
    stop(paste0(m.name, "  must be have ", nrow, " rows."))
  }
  if(ncol(m) != ncol){
    stop(paste0(m.name, "  must be have ", ncol, " columns"))
  }
  if (pos == TRUE){
    if(any(m < 0)){
      stop(paste0("Each element of ", m.name, " must be greater than or equal to 0."))
    }
  }
}

#' Check lifetable
#'
#' Check argument that passes a lifetable is correct. 
#' 
#' 
#' @return Error message if incorrect arguments are passed.
#'  
#' @keywords internal
check_lifetable <- function(lt){
  lt.name <- deparse(substitute(lt))
  if(is.null(lt$age)){
    stop(paste0("The age column in ", lt.name, " is missing."))
  }
  if(is.null(lt$qx)){
    stop(paste0("The qx column in ", lt.name, " is missing."))
  }
  if(lt$age[1] != 0){
    stop(paste0("Age in ", lt.name, " must start at 0."))
  }
  if(any(diff(lt$age) !=1) | !is.integer(lt$age)){
    stop(paste0(lt.name, " must report mortality rates by single-year of age."))
  }
  if(nrow(lt) != lt$age[length(lt$age)] + 1){
    stop("Each row in ", lt.name, " must represent a single year of age.")
  }
}

#' Check time to discontinuation
#'
#' Check that time to discontinuation argument in \link{sample_pars} is correct. 
#' 
#' 
#' @return Error message if incorrect arguments are passed.
#'  
#' @keywords internal
check_ttd <- function(pars, x_ttd){
  ttd.name <- deparse(substitute(pars))
  x.ttd.name <- deparse(substitute(x_ttd))
  names.dist <- c("exponential", "exp", "weibull", "gompertz", "gamma", "llogis",
                  "lnorm", "gengamma")
  if (!is.list(pars)){
      stop(paste0(ttd.name, " must be a list"))
  } 
  if (!all(sapply(pars, is.list))){
      stop(paste0("Some distributions in ", ttd.name, " are not lists"))
  }
  if (any(!names(pars) %in% names.dist)){
      stop(paste0("Distributions allowed in ", ttd.name, " are exponential, exp, ",
                "weibull, gompertz, gamma, llogis, lnorm, and gengamma."))
  }
  for (i in 1:length(names(pars))){ # loop over distributions
    dist <- names(pars)[i]
    pars.dist <- pars[[dist]]
    if (is.null(pars.dist$est)){
        stop(paste0("Element est is missing from the distribution ", dist, " in ", ttd.name))
    }
    if (is.null(pars.dist$vcov)){
        stop(paste0("Element vcov is missing from the distribution ", dist, " in ", 
                    ttd.name))
    }
    if (is.null(pars.dist$loc.index)){
       stop(paste0("Element loc.index is missing from the distribution ", dist, " in ",
                   ttd.name))
    } 
    if (is.null(pars.dist$anc1.index)){
      stop(paste0("Element anc1.index is missing from the distribution ", dist, " in ",
                  ttd.name))
    }
    if (is.null(pars.dist$anc2.index)){
        stop(paste0("Element anc2.index is missing from the distribution ", dist, " in ",
                    ttd.name))
    }
      
    if (dist %in% c("exponential", "exp")){
        if(!is.na(pars.dist$anc1.index)){
          stop(paste0("Element anc1.index from the distribution ", dist, " in ", ttd.name,
                     " must be NA."))
        }
        if(!is.na(pars.dist$anc2.index)){
          stop(paste0("Element anc2.index from the distribution ", dist, " in ", ttd.name,
                      " must be NA."))
        }
    }
    if (dist != "gengamma"){
      if(!is.na(pars.dist$anc2.index)){
        stop(paste0("Element anc2.index from the distribution ", dist, " in ", ttd.name,
                    " must be NA."))
      }
    } 
   check_vector(pars.dist$est)
   check_matrix(pars.dist$vcov, nrow = length(pars.dist$est), ncol = length(pars.dist$est))
   if (ncol(x_ttd) != length(pars.dist$loc.index)){
      stop(paste0("Number of columns in ", x.ttd.name, " is not equal to the length of ", 
                 "loc.index in ", ttd.name, " when using the ", dist, " distribution."))
   }
  } # end loop over distributions
}

#' Check NMA effect reduction
#'
#' Check argument that are passed to \link{sample_pars} related to the reduction
#' in treatment effects for bDMARD experienced patients.
#' 
#' 
#' @return Error message if incorrect arguments are passed.
#'  
#' @keywords internal
check_nma_k <- function(lower, upper){
  lower.name <- deparse(substitute(lower))
  upper.name <- deparse(substitute(upper))
  check_scalar(lower)
  check_scalar(upper)
  if (lower > upper){
    stop(paste0(upper.name, " must be greater than or equal to ", lower.name))
  }
}

