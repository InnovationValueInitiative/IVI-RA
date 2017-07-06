#' Sample patients
#'
#' Sample patients for rheumatoid arthritis individual patient simulation.
#' 
#' @param n Number of samples.
#' @param type Should male and female patients be heterogeneous or homogeneous. Default is homogeneous.
#' @param age_mean Mean age.
#' @param age_sd Standard deviation of age.
#' @param male_prop Proportion male.
#' @param bhaq0_mean Mean baseline HAQ score.
#' @param bhaq0_sd Standard deviation of baseline HAQ score.
#' @param wtmale Male weight.
#' @param wtfemale Female weight.
#' @param dis_dur_mean Mean disease duration.
#' @param dis_dur_sd Standard deviation of disease duration.
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
#' @return Matrix of patient characteristics. One row for each patient and one column
#' for each variable. Current variables are:
#' \describe{
#'   \item{age}{age in years}
#'   \item{male}{1 = male, 0 = female}
#'   \item{haq0}{Baseline HAQ score}
#'   \item{weight}{Patient weight}
#' }
#' @export
sample_pats <- function(n = 1, type = c("homog", "heterog"), age_mean = 55, age_sd = 13, male_prop = .21,
                      bhaq0_mean = 1.5, bhaq0_sd = 0.7, wtmale = 89, wtfemale = 75,
                      dis_dur_mean = 18.65, dis_dur_sd = 12.25, 
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
      haq0 <- rep(bhaq0_mean, n)
      dis_dur <- rep(dis_dur_mean, n)
      prev_dmards <- rep(prev_dmards_mean, n)
      das28 <- rep(das28_mean, n)
      sdai <- rep(sdai_mean, n)
      cdai <- rep(cdai_mean, n)
      da <- matrix(c(das28, sdai, cdai, haq0), ncol = 4)
  } else if (type == "heterog"){
      age <- msm::rtnorm(n, age_mean, age_sd, lower = 18, upper = 85) 
      haq0 <- msm::rtnorm(n, bhaq0_mean, bhaq0_sd, lower = 0, upper = 3)
      dis_dur <- msm::rtnorm(n, dis_dur_mean, dis_dur_sd, lower = 0)
      prev_dmards <- round(msm::rtnorm(n, prev_dmards_mean, prev_dmards_sd, lower = 0), 0)
      covmat <- matrix(0, nrow = 4, ncol = 4)
      diag(covmat) <- c(das28_sd^2, sdai_sd^2, cdai_sd^2, bhaq0_sd^2)
      covmat[1, 2] <- covmat[2, 1] <- cor_das28_sdai * das28_sd * sdai_sd
      covmat[1, 3] <- covmat[3, 1] <- cor_das28_cdai * das28_sd * cdai_sd
      covmat[1, 4] <- covmat[4, 1] <- cor_das28_haq * das28_sd * bhaq0_sd
      covmat[2, 3] <- covmat[3, 2] <- cor_sdai_cdai * sdai_sd * cdai_sd
      covmat[2, 4] <- covmat[4, 2] <- cor_sdai_haq * sdai_sd * bhaq0_sd
      covmat[3, 4] <- covmat[4, 3] <- cor_cdai_haq * cdai_sd * bhaq0_sd
      mu <- c(das28_mean, sdai_mean, cdai_mean, bhaq0_mean)
      da <- tmvtnorm::rtmvnorm(n, mean = mu, sigma = covmat,
                               lower = rep(0, length(mu)),
                               upper = c(9.4, 86, 76, 3))
  }
  x <- matrix(c(age, male, weight, dis_dur, prev_dmards, 
                da[, 1], da[, 2], da[, 3], da[, 4]), 
              nrow = n, ncol = 9, byrow = F)
  colnames(x) <- c("age", "male", "weight", "dis_dur", "prev_dmards", 
                   "das28", "sdai", "cdai", "haq0")
  return(x)
}

#' Lifetable for simulation
#' 
#' Generate lifetable matrix for use in simulation
#' 
#' @param ltmale
#' @param ltfemale
#' 
#' @return A list containing lifetables for females and males with qx transformed using the logit function. 
#' The lifetable has three columns: age, qx, and logit_qx.
#' 
#' @export
lt_data <- function(ltfemale, ltmale){
  ltmale <- ltmale[, .(age, qx)]
  ltfemale <- ltfemale[, .(age, qx)]
  ltmale[, logit_qx := ifelse(qx %in% c(0, 1), NA, log(qx/(1-qx)))]
  ltfemale[, logit_qx := ifelse(qx %in% c(0, 1), NA, log(qx/(1-qx)))]
  return(list(female = as.matrix(ltfemale), male = as.matrix(ltmale)))
}

#' Input data for simulation
#' 
#' Generate data inputs for the IPS.
#' 
#' @param patdata Matrix of patient data. Must contain variables generated from \link{sample_pats}:
#'  'age' for age, 'haq0' for baseline HAQ, 'male' as a indicator equal to
#' 1 if the patient is male and 0 if female, 'weight' for patient weight, 'dis_dur' for disease
#' duration and 'prev_dmards' for number of previous DMARDs. 
#' @param vars_mort A chararacter vector of variables in \code{patdata} to be used to adjust
#'  'qx' (using odds ratios) in the lifetables. Use '1' to include an intercept. 
#' @param vars_dur A character vector of variables in \code{patdata} to be used to predict treatment
#'  duration. Use '1' to include an intercept. 
#' 
#' @return A list containing the following data inputs:
#' \describe{
#'   \item{n}{Number of simulated patients}
#'   \item{haq0}{A vector of patient HAQ at baseline.}
#'   \item{age}{A vector of patient age at baseline.}
#'   \item{male}{A vector of patient gender (1 = male, 0 = female).}
#'   \item{dis.dur}{A vector of disease duration}
#'   \item{prev.dmards}{A vector of the number of previous DMARDs}
#'   \item{x.mort}{Design matrix for mortality adjustment with odds ratios}
#'   \item{x.dur}{Design matrix for treatment duration.}
#' }
#' 
#' @export
input_data <- function(patdata, vars_mort = "haq0", vars_dur = "1"){
  if ("1" %in% c(vars_mort, vars_dur)){
    int <- rep(1, nrow(patdata))
    patdata <- cbind(int, patdata)
    colnames(patdata)[1] <- "1"
  }
  x.mort <- patdata[, vars_mort, drop = FALSE]
  x.dur <- patdata[, vars_dur, drop = FALSE]
  return(list(n = nrow(patdata), haq0 = patdata[, "haq0"], age = patdata[, "age"],
              male = patdata[, "male"], weight = patdata[, "weight"],
              dis.dur = patdata[, "dis_dur"], prev.dmards = patdata[, "prev_dmards"],
              x.mort = x.mort, x.dur = x.dur))
}

