#' Sample a patient population
#'
#' Sample a patient population for use in the individual patient simulation (\link{sim_iviRA}).
#' 
#' @param n Number of samples.
#' @param type Should male and female patients be heterogeneous or homogeneous.
#' Default is homogeneous.
#' @param age_mean Mean age.
#' @param age_sd Standard deviation of age.
#' @param male_prop Proportion male.
#' @param haq0_mean Mean baseline HAQ (i.e., HAQ at the start of the model) score.
#' @param haq0_sd Standard deviation of baseline HAQ score.
#' @param wtmale Male weight.
#' @param wtfemale Female weight.
#' @param prev_dmards_mean Mean number of previous DMARDs.
#' @param prev_dmards_sd Standard deviation of number of previous DMARDs.
#' @param das28_mean Mean of DAS28.
#' @param das28_sd Standard deviation of DAS28.
#' @param sdai_mean Mean of SDAI.
#' @param sdai_sd Standard deviation of SDAI.
#' @param cdai_mean Mean of CDAI.
#' @param cdai_sd Standard deviation of CDAI.
#' @param cor_das28_sdai Correlation between DAS28 and SDAI.
#' @param cor_das28_cdai Correlation between DAS28 and CDAI.
#' @param cor_das28_haq Correlation between DAS28 and baseline HAQ.
#' @param cor_sdai_cdai Correlation between SDAI and CDAI.
#' @param cor_sdai_haq Correlation between SDAI and HAQ.
#' @param cor_cdai_haq Correlation between CDAI and HAQ.
#' 
#' 
#' @return Matrix of patient characteristics. One row for each patient and one column
#' for each variable. Current variables are:
#' \describe{
#'   \item{age}{Age in years.}
#'   \item{male}{1 = male, 0 = female.}
#'   \item{weight}{Patient weight in KG.}
#'   \item{prev_dmards}{Number of previous DMARDs.}
#'   \item{das28}{DAS28 score.}
#'   \item{sdai}{SDAI score.}
#'   \item{cdai}{CDAI score.}
#'   \item{haq0}{Baseline HAQ score.}
#' }
#' 
#' @examples
#' sample_pop(n = 10, type = "heterog", age_mean = 50)
#' 
#' @export
sample_pop <- function(n = 1, type = c("homog", "heterog"), age_mean = 55, age_sd = 13, male_prop = .21,
                      haq0_mean = 1.5, haq0_sd = 0.7, wtmale = 89, wtfemale = 75,
                      prev_dmards_mean = 3.28, prev_dmards_sd = 1.72,
                      das28_mean = 6, das28_sd = 1.2, 
                      sdai_mean = 43, sdai_sd = 13,
                      cdai_mean = 41, cdai_sd = 13,
                      cor_das28_sdai = .86, cor_das28_cdai = .86, cor_das28_haq = .38,
                      cor_sdai_cdai = .94, cor_sdai_haq = .34,
                      cor_cdai_haq = .34){
  if (age_mean > 85 | age_mean < 18){
    stop("Patients age must be between 18 and 85")
  }
  type <- match.arg(type)
  male <- rbinom(n, 1, male_prop)
  weight <- ifelse(male == 1, wtmale, wtfemale)
  if (type == "homog"){
      age <- rep(age_mean, n)
      haq0 <- rep(haq0_mean, n)
      prev_dmards <- rep(round(prev_dmards_mean, 0), n)
      das28 <- rep(das28_mean, n)
      sdai <- rep(sdai_mean, n)
      cdai <- rep(cdai_mean, n)
      da <- matrix(c(das28, sdai, cdai, haq0), ncol = 4)
  } else if (type == "heterog"){
      age <- msm::rtnorm(n, age_mean, age_sd, lower = 18, upper = 85) 
      haq0 <- msm::rtnorm(n, haq0_mean, haq0_sd, lower = 0, upper = 3)
      prev_dmards <- round(msm::rtnorm(n, prev_dmards_mean, prev_dmards_sd, lower = 0), 0)
      covmat <- matrix(0, nrow = 4, ncol = 4)
      diag(covmat) <- c(das28_sd^2, sdai_sd^2, cdai_sd^2, haq0_sd^2)
      covmat[1, 2] <- covmat[2, 1] <- cor_das28_sdai * das28_sd * sdai_sd
      covmat[1, 3] <- covmat[3, 1] <- cor_das28_cdai * das28_sd * cdai_sd
      covmat[1, 4] <- covmat[4, 1] <- cor_das28_haq * das28_sd * haq0_sd
      covmat[2, 3] <- covmat[3, 2] <- cor_sdai_cdai * sdai_sd * cdai_sd
      covmat[2, 4] <- covmat[4, 2] <- cor_sdai_haq * sdai_sd * haq0_sd
      covmat[3, 4] <- covmat[4, 3] <- cor_cdai_haq * cdai_sd * haq0_sd
      mu <- c(das28_mean, sdai_mean, cdai_mean, haq0_mean)
      da <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = covmat,
                               lower = rep(0, length(mu)),
                               upper = c(9.4, 86, 76, 3))
  }
  x <- matrix(c(age, male, weight, prev_dmards, 
                da[, 1], da[, 2], da[, 3], da[, 4]),
              nrow = n, ncol = 8, byrow = F)
  colnames(x) <- c("age", "male", "weight", "prev_dmards", 
                   "das28", "sdai", "cdai", "haq0")
  return(x)
}

#' Sample disease activity levels and HAQ 
#'
#' Sample disease activity levels and HAQ from a truncated multivariate normal distribution.
#' 
#' @param n Number of samples.
#' @param mean Mean vector.
#' @param sigma Covariance matrix.
#' @param lower Vector of lower truncation points.
#' @param upper Vector of upper truncation points. 
#' @keywords internal
sample_rtmvnorm_da <- function(n, mean, sigma, lower, upper){
  da <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = covmat,
                           lower = rep(-1, length(mu)),
                           upper = c(9.4, 86, 76, 3),
                           algorithm = "rejection")
}

#' Input data for simulation
#' 
#' Generate data inputs for the the individual patient simulation (\link{sim_iviRA}).
#' 
#' @param pop The patient population. A matrix that must contain variables generated from 
#' \link{sample_pop}: age' for age, 'haq0' for baseline HAQ, 'male' as a indicator equal to
#' 1 if the patient is male and 0 if female, 'weight' for patient weight, 'prev_dmards' 
#' for number of previous DMARDs, 'das28' for the patient's DAS28 score, 'sdai' for the 
#' patient's SDAI score, and 'cdai' for the patient's CDAI score. 
#' @param x_acr Design matrix where each column is a variable known at baseline that is 
#'   used to predict the relative treatment effects for ACR response. By default, 
#'   only includes an intercept, which implies  that there are no treatment-by-covariate 
#'   interactions. 
#' @param x_haq Design matrix where each column is a variable known at baseline that is
#'    used to predict the relative treatment effects for change in HAQ at 6 months from baseline. 
#'    By default, only includes an intercept, which implies that there are no 
#'    treatment-by-covariate interactions.
#' @param x_das28 Design matrix where each column is a variables known at baseline that 
#'   is used to predict the relative treatment effects for change in DAS28 at 6 months from 
#'   baseline. By default, only includes an intercept, which implies that there are no
#'    treatment-by-covariate interactions.
#' @param x_ttd_all Design matrix where each column is a variable influencing treatment duration
#'    in a model representative of all patients. The impact of each variable is determined
#'    by the sampled values of the coefficients used to predict the location parameter in 
#'    \code{ttd.all} returned by \link{sample_pars}.
#' @param x_ttd_da Design matrix where the first column is the intercept, the second column
#'   is a dummy variable used to indicate whether a patient has moderate disease activity, and
#'   the third column is dummy variable used to indicate whether a patient has high disease 
#'   activity (note, however, that the second and third columns are updated during the simulation
#'   according to the simulated disease activity level). All remaining columns after the 3rd column
#'   are variables influencing treatment duration. The impact of each variable is determined 
#'   by the sampled values of the coefficients used to predict the location parameter in 
#'    \code{ttd.da} returned by \link{sample_pars}.
#' @param x_ttd_eular Design matrix where each column is a variable influencing
#'   treatment duration in models stratified by EULAR response. The impact of each variable
#'   is determined by the sampled values of the coefficients used to predict the location
#'    parameter in \code{ttd.eular} returned by \link{sample_pars}.
#' @param x_mort Design matrix where each column is a variable influencing mortality. The impact 
#'   of each variable is determined by the parameter vector \code{mort.logor} returned by
#'   \link{sample_pars}.
#' @param x_attr Design matrix where each column is a variable related to treatment attributes
#'   related to the processes of care influencing utility. The impact 
#'   of each variable is determined by the parameter vector \code{tx.attr.utility} returned by
#'   \link{sample_pars}
#' 
#' @details If a design matrix is set to NULL, then a single column of ones is returned. 
#' In other words, if a design matrix is not specified, then it is assumed that an intercept
#'  only model will be used. 
#' 
#' @return A list containing the following data inputs:
#' \describe{
#'   \item{n}{Number of simulated patients.}
#'   \item{haq0}{A vector of patient HAQ at baseline.}
#'   \item{age}{A vector of patient age at baseline.}
#'   \item{male}{A vector of patient gender (1 = male, 0 = female).}
#'   \item{prev.dmards}{A vector of the number of previous DMARDs.}
#'   \item{x.acr}{Equivalent to \code{x.acr} passed as an argument to the function.}
#'   \item{x.haq}{Equivalent to \code{x.haq} passed as an argument to the function.}
#'   \item{x.das28}{Equivalent to \code{x.das28} passed as an argument to the function.}
#'   \item{x.ttd.all}{Equivalent to \code{x.ttd.all} passed as an argument to the function.}
#'   \item{x.ttd.da}{Equivalent to \code{x.ttd.da} passed as an argument to the function.}
#'  \item{x.ttd.eular}{Equivalent to \code{x.ttd.eular} passed as an argument to the function.}
#'   \item{x.mort}{Equivalent to \code{x.mort} passed as an argument to the function.}
#'   \item{x.attr}{Equivalent to \code{x.attr} passed as an argument to the function.}
#' }
#' 
#' @examples 
#' pop <- sample_pop(n = 100)
#' input.dat <- get_input_data(pop)
#' names(input.dat)
#' head(input.dat$haq0)
#' head(input.dat$x.haq)
#' @export
get_input_data <- function(pop, 
                           x_acr = NULL, x_haq = NULL, x_das28 = NULL,
                           x_ttd_all = NULL, x_ttd_da = NULL, x_ttd_eular = NULL, 
                           x_mort = NULL, 
                           x_attr = iviRA::utility.tx.attr$x){
  npats <- nrow(pop)
  
  # check pop
  if (!is.matrix(pop)){
    stop("pop must be a matrix.")
  }
  if (!"age" %in% colnames(pop)){
    stop("The variable age is missing from pop.")
  }
  if (!"haq0" %in% colnames(pop)){
    stop("The variable haq0 is missing from pop.")
  }
  if (!"male" %in% colnames(pop)){
    stop("The variable male is missing from pop.")
  }
  if (!"weight" %in% colnames(pop)){
    stop("The variable weight is missing from pop.")
  }
  if (!"prev_dmards" %in% colnames(pop)){
    stop("The variable prev_dmards is missing from pop.")
  }
  if (!"das28" %in% colnames(pop)){
    stop("The variable das28 is missing from pop.")
  }
  if (!"sdai" %in% colnames(pop)){
    stop("The variable sdai is missing from pop.")
  }
  if (!"cdai" %in% colnames(pop)){
    stop("The variable cdai is missing from pop.")
  }
  stopifnot(is.matrix(x_acr) | is.data.frame(x_acr) | is.null(x_acr), 
            is.matrix(x_haq) | is.data.frame(x_haq) | is.null(x_haq),
            is.matrix(x_das28) | is.data.frame(x_das28) | is.null(x_das28),
            is.matrix(x_ttd_all) | is.data.frame(x_ttd_all) | is.null(x_ttd_all),
            is.matrix(x_ttd_da) | is.data.frame(x_ttd_da) | is.null(x_ttd_da),
            is.matrix(x_ttd_eular) | is.data.frame(x_ttd_eular) | is.null(x_ttd_eular),
            is.matrix(x_mort) | is.data.frame(x_mort) | is.null(x_mort),
            is.matrix(x_attr) | is.data.frame(x_attr))
  
  # mortality
  if (is.null(x_mort)){
      x.mort <- pop[, "haq0", drop = FALSE]
  } else{
      x_mort <- as.matrix(x_mort)
      if (nrow(x_mort) != npats){
          stop("Number of rows in 'x_mort' must equal number of simulated patients.")
      }
    x.mort <- x_mort
  }
  
  # ACR response matrix for treatment by covariate interactions for d's
  if (is.null(x_acr)){
    x.acr <- matrix(1, nrow = nrow(pop), ncol = 1)
  } else{
    x_acr <- as.matrix(x_acr)
    if (nrow(x_acr) != npats){
      stop("Number of rows in 'x_acr' must equal number of simulated patients.")
    }
    x.acr <- x_acr
  }
  
  # Change in HAQ matrix for treatment by covariate interactions for d's
  if (is.null(x_haq)){
    x.haq <- matrix(1, nrow = nrow(pop), ncol = 1)
  } else{
    x_haq <- as.matrix(x_haq)
    if (nrow(x_haq) != npats){
      stop("Number of rows in 'x_acr' must equal number of simulated patients.")
    }
    x.haq <- x_haq
  }
  
  # DAS28 matrix for treatment by covariate interactions for d's
  if (is.null(x_das28)){
    x.das28 <- matrix(1, nrow = nrow(pop), ncol = 1)
  } else{
    x_das28 <- as.matrix(x_das28)
    if (nrow(x_das28) != npats){
      stop("Number of rows in 'x_das28' must equal number of simulated patients.")
    }
    x.da28 <- x.das28
  }
  
  # time to treatment discontinuation
    if(is.null(x_ttd_all)){
        x.ttd.all <- matrix(1, nrow = nrow(pop), ncol = 1)
    } else{
      x_ttd_all <- as.matrix(x_ttd_all)
      if (nrow(x_ttd_all) != npats){
        stop("Number of rows in 'x_ttd_all' must equal number of simulated patients.")
      }
      x.ttd.all <- x_ttd_all
    }
  
    if(is.null(x_ttd_eular)){
      x.ttd.eular <- matrix(1, nrow = nrow(pop), ncol = 1)
    } else{
      x_ttd_eular <- as.matrix(x_ttd_eular)
      if (nrow(x_ttd_eular) != npats){
        stop("Number of rows in 'x_ttd_eular' must equal number of simulated patients.")
      }
      x.ttd.eular <- x_ttd_eular
    }
  
    if(is.null(x_ttd_da)){
      x.ttd.da <- matrix(c(1, 0, 0), nrow = nrow(pop), ncol = 3, byrow = TRUE)
    } else{
      x_ttd_da <- as.matrix(x_ttd_da)
      if (nrow(x_ttd_da) != npats){
        stop("Number of rows in 'x_ttd_da' must equal number of simulated patients.")
      }
      x.ttd.da <- x_ttd_da
    }
  
  # treatment attributes
  x.attr <- x_attr

  
  # combine
  l <- list(n = nrow(pop), haq0 = pop[, "haq0"], age = pop[, "age"],
              male = pop[, "male"], das28 = pop[, "das28"],
              sdai = pop[, "sdai"], cdai = pop[, "cdai"],
              weight = pop[, "weight"], prev.dmards = pop[, "prev_dmards"],
              x.mort = x.mort, x.acr = x.acr, x.haq = x.haq, x.das28 = x.das28,
              x.ttd.all = x.ttd.all, x.ttd.eular = x.ttd.eular, x.ttd.da = x.ttd.da,
            x.attr = x.attr)
  class(l) <- "input_data"
  return(l)
}

