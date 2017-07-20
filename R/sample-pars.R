#' Sample parameters
#'
#' Randomly draw parameters from their (joint) probability distribution.
#' 
#' @param n Size of the posterior sample.
#' @param rebound_lower The rebound is the increase in HAQ following treatment. It is defined as 
#' a proportion \emph{f} times the size of the inititial treatment response. \code{rebound_lower} defines the
#' lower bound for \emph{f}. Default is 0.7, which implies
#' that the rebound post treatment is 0.7 times the initial treatment effect.
#' @param rebound_upper \code{rebound_upper} defines the
#' upper bound for \emph{f}. Default is 1, which implies
#' that the rebound post treatment is the same as the initial treatment effect.
#' @param ltfemale Lifetable for women. Must contain column 'age' for single-year of age and 'qx' for
#' the probability of death at a given age. 
#' @param ltmale Identical to \code{ltfemale} but for men.
#' @param acr2eular_mat A two-way frequency matrix with columns denoting EULAR response
#'  (none, moderate, good) and rows denoting ACR response  (<20, 20-50, 50-70, 70+).
#' @param treat_cost Treatment cost matrix in format of therapy.pars$cost.
#' @param mort_logor Log odds ratio of impact of baseline HAQ on probability of mortality.
#' @param mort_logor_se Standard error of log odds ratio of impact of baseline HAQ on probability of mortality.
#' @param mort_loghr_haqdif Log hazard ratio of impact of change in HAQ from baseline on mortality rate. A vector with
#' each element denoting (in order) hazard ratio for months 0-6, >6 - 12, >12 - 24, >24 -36, >36.
#' @param mort_loghr_se_haqdif Standard error of log hazard ratio of impact of change in HAQ from baseline on mortality rate.
#' @param ttd_all A list containing treatment duration parameters representative of all patients (i.e., unstratified).
#'  See 'Treatment duration'.
#' @param ttd_da A list containing treatment duration parameters. Covariates for moderate and high
#' disease activity. See 'Treatment duration'.
#' @param ttd_eular A list containing treatment duration parameters stratified by EULAR response. See 'Treatment duration'.
#' @param nma_acr_mean Posterior means for ACR response NMA parameters on probit scale for 
#' biologic naive patients (i.e., 1st line). ACR response is modeled using an ordered probit model.
#' @param nma_acr_vcov Variance-covariance matrix for ACR response NMA parameters on probit scale for 
#' biologic naive patients (i.e., 1st line). ACR response is modeled using an ordered probit model.
#' @param nma_acr_rr_lower Lower bound for proportion reduction (i.e., relative risk) in 
#' overlapping ACR response probabilities (ACR20/50/70) for biologic experienced patients 
#' (i.e., lines 2 and later).
#' @param nma_acr_rr_upper Upper bound for proportion reduction (i.e., relative risk) in
#'  overlapping ACR response probabilities (ACR20/50/70) for biologic experienced patients 
#'  (i.e., lines 2 and later).
#' @param nma_das28_mean Posterior means for DAS28 NMA parameters for biologic naive 
#' patients (i.e., 1st line). Change in DAS28 from baseline is modeled using a linear model.
#' @param nma_das28_vcov Variance-covariance matrix for DAS28 NMA paramters for biologic naive
#' patients (i.e., 1st line). Change in DAS28 from baseline is modeled using a linear model.
#' @param nma_das28_rr_lower Lower bound for proportion reduction (i.e., relative risk) in DAS28
#' for biologic experienced patients (i.e., lines 2 and later).
#' @param nma_das28_rr_upper Upper bound for proportion reduction (i.e., relative risk) in DAS28
#' for biologic experienced patients (i.e., lines 2 and later).
#' @param nma_haq_mean Posterior means for HAQ NMA parameters for biologic naive 
#' patients (i.e., 1st line). Change in HAQ from baseline is modeled using a linear model.
#' @param nma_haq_vcov Variance-covariance matrix for HAQ NMA paramters for biologic naive
#' patients (i.e., 1st line). Change in HAQ from baseline is modeled using a linear model.
#' @param nma_haq_rr_lower Lower bound for proportion reduction (i.e., relative risk) in HAQ
#' for biologic experienced patients (i.e., lines 2 and later).
#' @param nma_haq_rr_upper Upper bound for proportion reduction (i.e., relative risk) in HAQ
#' for biologic experienced patients (i.e., lines 2 and later).
#' @param treat_hist Is patient biologic naive or experienced. If naive NMA results for biologic 
#' naive patients are used during 1st line and biologic experienced NMA results are used
#' for subsequent lines. If experienced, NMA results for biologic
#' experienced patients are using during 1st and subsequent lines.
#' @param haq_lprog_therapy_mean Point estimate of linear yearly HAQ progression rate by therapy.
#' @param haq_lprog_therapy_se Standard error of linear yearly HAQ progression rate by therapy.
#' @param eular2haq_mean Mean HAQ change by Eular response category.
#' @param eular2haq_se Standard error of mean HAQ change by Eular response category.
#' @param acr2haq_mean Mean HAQ change by ACR response category.
#' @param acr2haq_se Standard error of mean HAQ change by ACR response category.
#' @param acr2sdai_lower Lower bound for change in SDAI by ACR response category.
#' @param acr2sdai_upper Upper bound for change in SDAI by ACR response category.
#' @param acr2cdai_lower Lower bound for change in CDAI by ACR response category.
#' @param acr2cdai_upper Upper bound for change in CDAI by ACR response category.
#' @param acr2das28_lower Lower bound for change in DAS28 by ACR response category.
#' @param acr2das28_upper Upper bound for change in DAS28 by ACR response category.
#' @param haq_lprog_age_mean Impact of age on annual linear haq progression rate.
#' @param haq_lprog_age_se Standard error of impact of age on annual linear haq progression rate.
#' @param hosp_days_mean Vector denoting average number of hospital days for HAQ < 0.5,
#' 0.5 <= HAQ < 1, 1 <= HAQ < 1.5, 1.5 <= HAQ < 2, 2 <= HAQ < 2.5, HAQ >= 2.5. 
#' @param hosp_days_se Vector denoting standard error of average number of hospital days for HAQ < 0.5,
#' 0.5 <= HAQ < 1, 1 <= HAQ < 1.5, 1.5 <= HAQ < 2, 2 <= HAQ < 2.5, HAQ >= 2.5. 
#' @param hosp_cost_mean Mean of daily hospital cost.
#' @param hosp_cost_se Standard error of dail hospital cost.
#' @param mgmt_cost_mean Mean of costs of services (in order: chest x-ray, x-ray visit,
#'  outpatient followup, Mantoux tuberculin skin test) for general management of RA.
#' @param mgmt_cost_se Standard error of mean of costs of services for general management of RA
#'  for general management of RA.
#' @param pl_mean Mean annual productivity loss per 1-unit increase in HAQ.
#' @param pl_se Standard error of mean annual productivity loss per 1-unit increase in HAQ.
#' @param ttsi Paramters of survival model used to estimate time to serious infection.
#' @param si_ul One month loss in utility from a serious infection. 
#' @param si_ul_range Range used to vary serious infection utility loss. Default is to calculate upper and lower bound by multiplying 
#' \code{si_ul} by 1 +/- 0.2 (i.e. a 20\% change).
#' @param si_cost Cost of a serious infection.
#' @param si_cost_range Range used to vary serious infection cost. Default is to calculate upper and lower bound by multiplying 
#' \code{si_cost} by 1 +/- 0.2 (i.e. a 20\% change).
#' @param therapy_names Vector of therapy names.
#' 
#' @return List containing samples for the following model parameters:
#' 
#' \describe{
#'  \item{treat.cost}{\code{data.table} of treatment cost data.}
#'  \item{ltmale}{Lifetable for males with qx transformed using the logit function.}
#'   \item{ltfemale}{Lifetable for females with qx transformed using the logit function.}
#'  \item{acr2eular}{An array of matrices. Each matrix represents a random sample of the conditional 
#'  probability of each EULAR response category for a given ACR response.}
#'   \item{logor.mort}{Matrix of log odds ratio used to adjust mortality. One row for each sample
#'   and one column for each variable used to adjust mortality.}
#'   \item{mort.loghr.haqdif}{Matrix of the log hazard ratio of the impact of a change in HAQ from baseline on mortality. Columns denote
#'   hazard ratios at times < 6 months, months 6 - <12, months 12 - <24, months 24 - <36, and months 36+.}
#'   \item{ttd.eular.mod}{Matrix of coefficients from a survival model. From a survival model for moderate Eular responders}
#'   \item{ttd.eular.good}{Matrix of coefficients for the location parameter and a vector
#'    of sampled values of the ancillary parameter. From a survival model for good Eular responders.}
#'   \item{acr}{A list containing six objects. \code{p1} is an array containing information on the
#'   probability of achieving the four mutually exclusive ACR categories (ACR<20, ACR 20-50,
#'   ACR 50-70, ACR70+). The array stores matrices with each row each row a sampled parameter value
#'   and columns denoting the mutuaully exclusive ACR categories. Therapies are indexed using
#'   the third dimension. \code{p1.overalap} is identical to \code{p1} but for overlapping ACR
#'   response categories (ACR <20, ACR 20+, ACR 50+, ACR 70+). \code{p1} and \code{p1.overlap} are
#'   for the first treatment line. \code{p2} and \code{p2.overlap} are the same as 
#'   \code{p1} and \code{p1.overlap} but for treatments second line and later. \code{rr} is the 
#'   sampled values of the probabilility of overalapping ACR response categories for biologic
#'   experienced patients relative to the probability for biologic naive patients (i.e., 
#'   the relative risk). Finally, \code{pars} is a matrix of sampled values of model coefficients 
#'   from the ordered logistic regression used for the NMA.} 
#'   \item{das28}{A list of four elements. \code{dy1} and \code{dy2} are matrices with each row 
#'   a sampled parameter value and columns denoting the mean change in DAS28 for each
#'   therapy. \code{dy1} is for first line treatment and \code{dy2} is for second line treatment
#'   and later. \code{pars} is a matrix of sampled values of model coefficients from the linear
#'  model used for the NMA. Finally, \code{rr} is the sampled values of the mean change in
#'    DAS28 from baseline for biologic experienced patients as a fraction of the mean change for
#'    biologic naive patients.}
#'   \item{haq}{Identical to DAS28 but for the HAQ score.}
#'   \item{eular2haq}{A matrix of sampled HAQ changes by Eular response category. The matrix has
#'    three columns for no response, moderate response, and good response.}
#'    \item{acr2haq}{A matrix of sampled HAQ changes by ACR response category. The matrix has
#'    four columns for ACR < 20, ACR 20-50, ACR 50-70, and ACR 70+.}
#'   \item{haq.lprog.therapy}{A matrix of sampled yearly linear change in HAQ by therapy. The matrix has one column
#'    for each therapy in \code{therapy.pars}.}
#'    \item{haq.lprog.age}{A matrix of sampled yearly linear change in HAQ by age. The matrix
#'    has three columns for age < 40, age 40-64, and age 65+.}
#'    \item{haq.lcgm}{A list of two elements containing parameters from the latent class growth
#'     model. The first element is \code{delta} which is a an array of sampled matrices with
#'      each matrix containing coefficients predicting class membership. Rows are classes and columns index
#'     variables. \code{beta} is similar to \code{delta}, but each matrix contains coefficients 
#'     predicting HAQ as a function of time using a quadratic polynomial model.
#'    }
#'    \item{mixture.utility}{A list containing samples of all parameters in the Hernandez Alva (2013) mixture model. See 'Sampled mixture model parameters'
#'    for details.}
#'    \item{wailoo.utility}{A matrix of sampled regression coefficients from the model mapping HAQ to EQ5D utility in Wailoo (2006). Variables are
#'    (in order) "int" (intercept), "age" (patient age), "dis_dur" (disease duration), "haq0" (baseline HAQ), "male" (1 = male, 0 = female),
#'    "prev_dmards" (number of previous DMARDs), and "haq" (current HAQ).}
#'    \item{hosp.cost}{A list of two matrices \code{hosp.days} and \code{cost.pday}. \code{hosp.days} is sample of hospital days by HAQ category; the 
#'    columns of the matrix are the six HAQ categories (HAQ < 0.5, 0.5 <= HAQ < 1, 1 <= HAQ < 1.5, 1.5 <= HAQ < 2, 2 <= HAQ < 2.5, HAQ >= 2.5). 
#'    in \code{hosp.days} are HAQ. \code{cost.pday} is a sample of the costs per hospital day by HAQ category; the columns are the same six HAQ
#'    categories as in \code{hosp.days}.}
#'    \item{mgmt.cost}{Matrix of sampled values of general management costs. Each column is a different category of costs (
#'    chest x-ray, x-ray visit, outpatient follow-up, and Mantoux tuberculin skin test). }
#'    \item{prod.loss}{Vector of sampled values of decrease in wages (e.g. productivity loss) per unit increase in HAQ.}
#'    \item{ttsi}{Matrix of coefficients from a time-to event model predicting time to serious infection.
#'     Columns numbers coincide with therapy indices.}
#'    \item{si.cost}{Vector of sampled values of the medical cost of a serious infection.}
#'    \item{si.ul}{Vector of the sampled values of the annualized utility loss from a serious infection.}
#' }
#'
#' @section Time to treatment discontinuation:
#' Time to treatment discontinuation parameters should be contained in a list of lists. The top-level
#'  list identifies the name of the probability distribution; the possible distributions are the 
#'  exponential (\code{exponential}), Weibull (\code{weibull}), Gompertz (\code{gompertz}),
#' gamma (\code{gamma}), log-logistic (\code{llogis}), lognormal (\code{lnorm}), and generalized gamma (\code{gengamma}). Each distribution 
#' should also contain a list with five elements:
#' \describe{
#' \item{est}{A vector of the maximum likelihood estimates of the parameters. }
#' \item{vcov}{The variance-covariance of the parameters.}
#' \item{loc.index}{A vector of the indices of the location parameters.}
#' \item{anc1.index}{A vector of the indices of the first ancillary parameter.}
#' \item{anc2.index}{A vector of indices of the second ancillary parameter.}
#' }
#' The maximum likelihood estimates should be transformed to the real line. For example, if the model is fit using 
#' \code{flexsurvreg} in the \code{flexsurv} package, the output should be returned from \code{res.t}.  
#'
#'@section Sampled mixture model parameters:
#' The sampled mixture model parameters are contained in a list containing the following:
#' \describe{
#' \item{beta1}{Coefficients for class 1 explanatory variables. A matrix of random draws where each
#'  column is an explanatory variable.}
#' \item{beta2}{Coefficients for class 2 explanatory variables. A matrix of random draws where each
#'  column is an explanatory variable.}
#' \item{beta3}{Coefficients for class 3 explanatory variables. A matrix of random draws where each
#'  column is an explanatory variable.}
#' \item{beta4}{Coefficients for class 4 explanatory variables. A matrix of random draws where each
#'  column is an explanatory variable.}
#' \item{alpha1}{Random effects intecept term for class 1. A vector of random draws.}
#' \item{alpha2}{Random effects intecept term for class 2. A vector of random draws.}
#' \item{alpha3}{Random effects intecept term for class 3. A vector of random draws.}
#' \item{alpha4}{Random effects intecept term for class 4. A vector of random draws.}
#' \item{alpha}{Random effects term for male indicator variable. A vector of random draws.}
#' \item{epsilon1}{Variance for class 1. A vector of random draws.}
#' \item{epsilon2}{Variance for class 2. A vector of random draws.}
#' \item{epsilon3}{Variance for class 3. A vector of random draws.}
#' \item{epsilon4}{Variance for class 4. A vector of random draws.}
#' \item{mu}{Random effects variance term. A vector of random draws.}
#' \item{delta}{Coefficients for explanatory variables explaining the probbaility
#' of class membership. An array of matrices where each matrix is a random draw. There are three rows in each matrix (one for each class) and
#' four columns (one for each explanatory variable).}
#' }
#' The list also contains summary statistics for pain and HAQ, which is needed to simulate pain. 
#' In particular, the object 'pain' in the list is a list containing:
#' \describe{
#' \item{pain.mean}{Mean of pain score in the population.}
#' \item{haq.mean}{Mean of HAQ score in the population.}
#' \item{pain.var}{Variance of pain score in the population.}
#' \item{haq.var}{Variance of HAQ in the population.}
#' \item{painhaq.cor}{Correlation between pain and HAQ in the population.}
#' }
#' 
#' @export
sample_pars <- function(n = 100, rebound_lower = .7, rebound_upper = 1,
                       ltfemale = iviRA::lifetable.female, ltmale = iviRA::lifetable.male,
                       acr2eular_mat = iviRA::acr2eular,
                       treat_cost = iviRA::therapy.pars$cost,
                       mort_logor = iviRA::mort.or$logor, mort_logor_se = iviRA::mort.or$logor_se,
                       mort_loghr_haqdif = iviRA::mort.hr.haqdif$loghr,
                       mort_loghr_se_haqdif = iviRA::mort.hr.haqdif$loghr_se,
                       ttd_all = iviRA::ttd.all, ttd_da = iviRA::ttd.da, ttd_eular = iviRA::ttd.eular,
                       nma_acr_mean = iviRA::therapy.pars$nma.acr.naive$mean,
                       nma_acr_vcov = iviRA::therapy.pars$nma.acr.naive$vcov,
                       nma_acr_rr_lower = .75, nma_acr_rr_upper = .92,
                       nma_das28_mean = iviRA::therapy.pars$nma.das28.naive$mean,
                       nma_das28_vcov = iviRA::therapy.pars$nma.das28.naive$vcov,
                       nma_das28_rr_lower = .75, nma_das28_rr_upper = .92,
                       nma_haq_mean = iviRA::therapy.pars$nma.haq.naive$mean,
                       nma_haq_vcov = iviRA::therapy.pars$nma.haq.naive$vcov,
                       nma_haq_rr_lower = .75, nma_haq_rr_upper = .92,
                       treat_hist = c("naive", "exp"),
                       haq_lprog_therapy_mean = iviRA::therapy.pars$haq.lprog$est,
                       haq_lprog_therapy_se = iviRA::therapy.pars$haq.lprog$se,
                       eular2haq_mean = iviRA::eular2haq$mean, 
                       eular2haq_se = iviRA::eular2haq$se,
                       acr2haq_mean = iviRA::acr2haq$mean,
                       acr2haq_se = iviRA::acr2haq$se,
                       acr2sdai_lower = iviRA::acr2sdai$inception$lower,
                       acr2sdai_upper = iviRA::acr2sdai$inception$upper,
                       acr2cdai_lower = iviRA::acr2cdai$inception$lower,
                       acr2cdai_upper = iviRA::acr2cdai$inception$upper,
                       acr2das28_lower = iviRA::acr2das28$inception$lower,
                       acr2das28_upper = iviRA::acr2das28$inception$upper,
                       haq_lprog_age_mean = iviRA::haq.lprog.age$est,
                       haq_lprog_age_se = iviRA::haq.lprog.age$se,
                       hosp_days_mean = c(.26, .13, .51, .72, 1.86, 4.16),
                       hosp_days_se = c(.5, .5, .5, .5, .5, .5),
                       hosp_cost_mean = rep(1251, 6),
                       hosp_cost_se = rep(191, 6),
                       mgmt_cost_mean = iviRA::mgmt.cost$est,
                       mgmt_cost_se = iviRA::mgmt.cost$se,
                       pl_mean = 5853.30, pl_se = 1526.691,
                       ttsi = iviRA::therapy.pars$si,
                       si_cost = 5873, si_cost_range = .2,
                       si_ul = .156, si_ul_range = .2,
                       therapy_names = iviRA::therapy.pars$info$sname,
                       util_mixture_pain = iviRA::pain){
  acr.cats <- c("ACR <20", "ACR 20-50", "ACR 50-70", "ACR 70+")
  sim <- list()
  sim$n <- n
  sim$rebound <- runif(n, rebound_lower, rebound_upper)
  sim$lt <- lt_data(ltmale, ltfemale)
  sim$treat.cost <- calc_treat_cost(treat_cost)
  sim$acr2eular <- sample_dirichlets(n, acr2eular_mat)
  sim$acr2sdai <- sample_uniforms(n, acr2sdai_lower, acr2sdai_upper, acr.cats)
  sim$acr2cdai <- sample_uniforms(n, acr2cdai_lower, acr2cdai_upper, acr.cats)
  sim$acr2das28 <- sample_uniforms(n, acr2das28_lower, acr2das28_upper, acr.cats)
  sim$logor.mort <- sample_mvnorm(n, mort_logor, mort_logor_se^2)
  sim$mort.loghr.haqdif <- sample_normals(n, mort_loghr_haqdif, mort_loghr_se_haqdif,
                                         col_names = paste0("month_", c("less6", "6to12", "12to24", "24to36", "36to48")))
  sim$ttd.all <- sample_survpars(n, ttd_all)
  sim$ttd.da <- sample_survpars(n, ttd_da)
  sim$ttd.eular <- sample_stratified_survpars(n, ttd_eular)
  treat_hist <- match.arg(treat_hist)
  sim$acr <- sample_nma_acr(n, nma_acr_mean, nma_acr_vcov, rr_lower = nma_acr_rr_lower,
                              rr_upper = nma_acr_rr_upper, hist = treat_hist)
  sim$das28 <- sample_nma_lm(n, nma_das28_mean, nma_das28_vcov, rr_lower = nma_das28_rr_lower,
                                 rr_upper = nma_das28_rr_upper, hist = treat_hist)
  sim$haq <- sample_nma_lm(n, nma_haq_mean, nma_haq_vcov, rr_lower = nma_haq_rr_lower,
                            rr_upper = nma_haq_rr_upper, hist = treat_hist)
  sim$eular2haq <- sample_normals(n, eular2haq_mean, eular2haq_se,
                                 col_names = c("no_response", "moderate_response", "good_response"))
  sim$acr2haq <- sample_normals(n, acr2haq_mean, acr2haq_se, acr.cats)
  sim$haq.lprog.therapy <- sample_normals(n, haq_lprog_therapy_mean,
                                                  haq_lprog_therapy_se)
  sim$haq.lprog.age <- sample_normals(n, haq_lprog_age_mean, haq_lprog_age_se,
                                     col_names =  c("age_less40", "age40to64", "age_65plus"))
  sim$haq.lcgm <- sample_pars_haq_lcgm(n)
  sim$mixture.utility <- sample_pars_utility_mixture(n, util_mixture_pain)
  sim$wailoo.utility <- sample_normals(n, iviRA::util.wailoo.pars$coef, iviRA::util.wailoo.pars$se,
                                      col_names = names(util.wailoo.pars$coef))
  hosp.cost.names <- c("haq_less0.5", "haq0.5to1", "haq1to1.5", "haq1.5to2", "haq2to2.5", "haq2.5plus")
  sim$hosp.cost <- list(hosp.days = sample_gammas(n, hosp_days_mean, hosp_days_se,
                                                 col_names = hosp.cost.names),
                     cost.pday = sample_gammas(n, mean = hosp_cost_mean, se = hosp_cost_se,
                                              col_names = hosp.cost.names))
  sim$mgmt.cost <- sample_gammas(n, mgmt_cost_mean, mgmt_cost_se,
                                col_names = c("chest_xray", "xray_visit", "outpatient_followup", "tuberculin_test"))
  sim$prod.loss <- rnorm(n, pl_mean, pl_se)
  sim$ttsi <- sample_survpars(n, ttsi)
  sim$si.cost <- runif(n, si_cost * (1 - si_cost_range),  si_cost * (1 + si_cost_range))
  sim$si.ul <- runif(n, si_ul * (1 - si_ul_range), si_ul * (1 + si_ul_range))
  return(sim)
}

#'  Calculate treatment cost
#'
#'  Calculate treatment cost from data of dosing and prices
#' 
#' @return Data table
#' 
#' @export
calc_treat_cost <- function(x){
  
  # cost first 6 months
  x[!sname %in% c("ifx", "abtsc"), ':=' (init_infusion_cost = init_num_doses * infusion_cost,
       init_rx_cost =  init_dose_val/strength_val * init_num_doses * wac_per_unit)]
  x[sname == "ifx", ':=' (init_infusion_cost = init_num_doses * infusion_cost, init_rx_cost = 0)]
  x.abtiv <- x[sname == "abtiv"]
  x[sname == "abtsc", ':=' (init_infusion_cost = 1 * x.abtiv$infusion_cost, 
                            init_rx_cost = 1 * (x.abtiv$init_dose_val/x.abtiv$strength_val * x.abtiv$wac_per_unit) + 
                              init_dose_val/strength_val * init_num_doses * wac_per_unit)]
  
  # annual cost 6+ months
  x[!sname %in% c("ifx"), ':=' (ann_infusion_cost = ann_num_doses * infusion_cost,
                                ann_rx_cost = ann_dose_val/strength_val * ann_num_doses * wac_per_unit)]
  x[sname == "ifx", ':=' (ann_infusion_cost = ann_num_doses * infusion_cost, ann_rx_cost = 0)]
  
  # placebo + non-biologic
  y <- x[, .(sname, init_infusion_cost, init_rx_cost, ann_infusion_cost, ann_rx_cost)]
  y <- rbind(y, data.table(sname = "placebo", init_infusion_cost = 0, init_rx_cost = 0, 
                           ann_infusion_cost = 0, ann_rx_cost = 0))
  y.cdmards <- y[sname == "cdmards"]
  y <- rbind(y, data.table(sname = "nbt", init_infusion_cost = y.cdmards$init_infusion_cost,
                           init_rx_cost = y.cdmards$init_rx_cost, 
                           ann_infusion_cost = y.cdmards$ann_infusion_cost, ann_rx_cost = y.cdmards$ann_rx_cost))
  
  # combination therapies
  bmtx <- c("abtiv", "ada", "etn", "gol", "ifx", "tcz", "czp", 
            "abtsc", "rtx", "tof")
  y2 <- data.table(sname = c("abtivmtx", "adamtx", "etnmtx", "golmtx", 
                             "ifxmtx", "tczmtx", "czpmtx", "abtscmtx", "rtxmtx", 
                             "tofmtx"))
  cnames <- c("init_infusion_cost", "init_rx_cost", "ann_infusion_cost", "ann_rx_cost")
  for (i in 1:length(bmtx)){
    for (j in cnames){
      val <- as.numeric(y[sname == bmtx[i], j, with = FALSE] + y[sname == "cdmards", j, with = F])
      y2[i, c(j) := val]
    }
  }
  y <- rbind(y, y2)
  y.tt <- y[sname %in% c("cdmards", "hcl", "ssz")]
  y <- rbind(y, data.table(sname = "tt", init_infusion_cost = sum(y.tt$init_infusion_cost),
                           init_rx_cost = sum(y.tt$init_rx_cost),
                           ann_infusion_cost = sum(y.tt$ann_infusion_cost),
                           ann_rx_cost = sum(y.tt$ann_rx_cost)))
  
  # subset/order therapies
  y <- y[sname %in% therapy.pars$info$sname]
  y <- merge(y, therapy.pars$info[, .(sname, mname)], by = "sname")
  y <- y[, .(sname, mname, init_infusion_cost, init_rx_cost, ann_infusion_cost, ann_rx_cost)]
  y <- y[order(match(sname, therapy.pars$info$sname))]
  
  # weight-based
  y[, weight_based := ifelse(sname == "ifx" | sname == "ifxmtx", 1, 0)]
  x.ifx <- x[sname == "ifx"]
  y[, init_wgt_slope := ifelse(sname == "ifxmtx", x.ifx$init_dose_val, 0)]
  y[, ann_wgt_slope := ifelse(sname == "ifxmtx", x.ifx$ann_dose_val, 0)]
  y[, init_util := ifelse(sname == "ifxmtx", x.ifx$init_num_doses, 0)]
  y[, ann_util := ifelse(sname == "ifxmtx", x.ifx$ann_num_doses, 0)]
  y[, strength := ifelse(sname == "ifxmtx", x.ifx$strength_val, 0)]
  y[, price := ifelse(sname == "ifxmtx", x.ifx$wac_per_unit, 0)]
  
  # return
  return(y)
}

#' Sample from beta distributions
#'
#' Sample from multiple independent beta distributions. Produces a unique sample of size \code{n} 
#' using each element in \code{shape1} and \code{shape2}. 
#' @param n Number of observations.
#' @param shape1 Alpha parameter in beta distribution.
#' @param shape2 Beta parameter in beta distribution.
#' @param col_names Column names to returned matrix.
#' 
#' @return Matrix with each column a sample from a beta distribution. 
#' 
#' @export
sample_betas <- function(n, shape1, shape2, col_names = NULL){
  s <- matrix(rbeta(n * length(shape1), shape1, shape2),
              nrow = n, ncol = length(shape1), byrow = T)
  colnames(s) <- col_names
  return(s)
}

#' Sample matrices using a Dirichlet distribution 
#'
#' A wrapper using the function \code{hesim::rdiricihlet_mat}. Maintains the row and column names of \code{mat}.
#' @param n Number of samples to draw.
#' @param mat A matrix where each row is a separte vector of shape parameters.
#' 
#' @return An array of matrices representing a sample from the dirichlet distribution for each row.
#' 
#' @export
sample_dirichlets <- function(n, mat){
    samp <- hesim::rdirichlet_mat(n, mat)
    dimnames(samp) <- list(rownames(mat),
                        colnames(mat),
                        seq(1, n))
  return(samp)
}

#' Sample from a multivariate normal distribution
#'
#' A wrapper for \code{MASS::mvrnorm} that returns a matrix (rather than a vector) if n = 1. Also,
#' elements in which the mean vector has missing values.   
#' 
#' @param n The number of samples required.
#' @param mu A vector giving the means of the variables.
#' @param Sigma A positive-definite symmetric matrix specifying the covariance matrix of the variables. 
#' 
#' @return Matrix of samples from the multivariate distribution.
#' 
#' @export
sample_mvnorm <- function(n, mu, Sigma){
  mu.na.pos <- which(is.na(mu))
  if (length(mu.na.pos) == 0){
      sim <- MASS::mvrnorm(n, mu, Sigma)
      if (class(sim) == "numeric") sim <- t(as.matrix(sim)) 
  } else{
      if (length(mu) == 1){
        stop("mu is missing and of length 1.")
      }
      mu.na.names <- names(mu)[mu.na.pos]
      mu.complete.pos <- which(!is.na(mu))
      mu <- mu[mu.complete.pos]
      Sigma <- Sigma[mu.complete.pos, mu.complete.pos]
      if (sum(is.na(Sigma)) > 0){
        stop("Variance-covariance matrix used to draw from the multivariate normal distribution
         contains missing values.")
      }
      sim <- MASS::mvrnorm(n, mu, Sigma)
      if (class(sim) == "numeric") sim <- t(as.matrix(sim)) 
      sim.na <- matrix(NA, nrow = n, ncol = length(mu.na.pos))
      colnames(sim.na) <- mu.na.names
      sim <- cbind(sim, sim.na)
      pos.sim <- c(mu.complete.pos, mu.na.pos)
      pos <- seq(1, ncol(sim))
      sim <- sim[, match(pos, pos.sim), drop = FALSE]
  }
  return(sim)
}

#' Sample from normal distributions 
#'
#' Sample from multiple independent normal distributions. Produces a unique sample of size \code{n} using each element in 
#' \code{mean} and \code{sd}. Equivalent to sampling from a multivariate normal distribution with no covariance. Each sample using \code{rnorm}. 
#' 
#' @param n Number of observations.
#' @param mean Vector of means.
#' @param sd Vector of stadnard deviations. 
#' @param col_names Column names to returned matrix.
#' 
#' @return Matrix with each column a sample from a normal distribution. 
#' 
#' @export
sample_normals <- function(n, mean, sd, col_names = NULL){
  s <- matrix(rnorm(n * length(mean), mean, sd),
              nrow = n, ncol = length(mean), byrow = T)
  colnames(s) <- col_names
  return(s)
}

#' Sample from uniform distributions 
#'
#' Sample from multiple independent uniform distributions. Produces a unique sample of size \code{n} using each element in 
#' \code{lower} and \code{upper}. Each sample using \code{runif}. 
#' 
#' @param n Number of observations.
#' @param lower Vector of lower bounds.
#' @param upper Vector of upper bounds. 
#' @param col_names Column names to returned matrix.
#' 
#' @return Matrix with each column a sample from a uniform distribution. 
#' 
#' @export
sample_uniforms <- function(n, lower, upper, col_names = NULL){
  s <- matrix(runif(n * length(lower), lower, upper),
              nrow = n, ncol = length(lower), byrow = T)
  colnames(s) <- col_names
  return(s)
}

#' Sample from Bayesian linear model
#'
#' Sample change in outcomes from Bayesian linear model for NMA.
#' 
#' @param nsims Number of observations.
#' @param m Mean for each coefficient.
#' @param vcov Variance-covariance matrix of coefficients.
#' @param rr_lower Lower bound for relative risk.
#' @param rr_upper Upper bound for relative risk.
#' @param hist Patient history equivalent to treat_hist in \link{sample_pars}
#' @return List containing posterior samples of changes in outcomes.
#' 
#' @export
sample_nma_lm <- function(nsims, m, vcov, rr_lower, rr_upper, hist){
  rr.sim <- runif(nsims, rr_lower, rr_upper)
  sim <- sample_mvnorm(nsims, m, vcov)
  if (hist == "naive"){
    dy1 <- nma_lm2prob(A = sim[, "A"], delta = sim[, 2:ncol(sim), drop = FALSE])
    dy2 <- dy1 * rr.sim
  } else if (hist == "exp"){
    dy1 <- dy2 <- nma_lm2prob(A = sim[, "A"], delta = sim[, 2:ncol(sim), drop = FALSE]) * rr.sim
  }
  colnames(dy1) <- colnames(dy2) <- therapy.pars$info$sname
  return(list(dy1 = dy1, dy2 = dy2, pars = sim, rr = rr.sim))
}

#' Sample ACR response from ordered probit NMA
#'
#' Sample ACR response from an NMA using a multinomial likelihood and a probit link.
#' See for example section 3.6 and example 6 in "Nice DSU Technical Support Document 2:
#' A Generalized Linear Modeling Framework for Pairwise and Network Meta-Analysis of
#' Randomized Controlled Trials"
#' 
#' @param nsims Number of observations.
#' @param m Mean for each coefficient.
#' @param vcov Variance-covariance matrix of coefficients.
#' @param rr_lower Lower bound for relative risk.
#' @param rr_upper Upper bound for relative risk.
#' @param hist Patient history equivalent to treat_hist in \link{sample_pars}.
#' @return List containing posterior sample of ACR response for each therapy
#' 
#' @export
sample_nma_acr <- function(nsims, m, vcov, rr_lower, rr_upper, hist){
  rr2.sim <- runif(nsims, rr_lower, rr_upper)
  if (hist == "naive"){
      rr1.sim <- 1
  } else if (hist == "exp"){
      rr1.sim <- rr2.sim
  }
  sim <- sample_mvnorm(nsims, m, vcov)
  p1 <- nma_acrprob(A = sim[, "A"], z2 = sim[, "z2"], z3 = sim[, "z3"],
                      delta = sim[, 5:ncol(sim), drop = FALSE],
                      rr = rr1.sim)
  p2 <- nma_acrprob(A = sim[, "A"], z2 = sim[, "z2"], z3 = sim[, "z3"],
                    delta = sim[, 5:ncol(sim), drop = FALSE],
                    rr = rr2.sim) 
  dimnames(p1$non.overlap)[[3]] <- dimnames(p2$non.overlap)[[3]] <- therapy.pars$info$sname
  dimnames(p1$overlap)[[3]] <- dimnames(p2$overlap)[[3]] <-  therapy.pars$info$sname
  return(list(p1 = p1$non.overlap, p1.overlap = p1$overlap,
              p2 = p2$non.overlap, p2.overlap = p2$overlap,
              rr = rr2.sim, pars = sim))
}

#' Sample survival parameters
#'
#' Generate a random sample of parameters for each survival distribution. Parameters are
#' sampled using a multivariate normal distribution.
#' 
#' @param nsims Size of the posterior sample.
#' @param x A list of survival parameters. For example, see 'Treatment duration' in \link{sample_pars}.
#' 
#' @return A list of matrices containing random draws of the parameters of the survival distribution. One matrix
#' for the location parameter and one matrix for each of the ancillary parameters. 
#' 
#' @export
sample_survpars <- function(nsims, x){
  dists <- names(x) 
  l <- vector(length(dists), mode = "list")
  names(l) <- dists
  for (i in 1:length(dists)){
    xi <- x[[dists[i]]]
    samp <- MASS::mvrnorm(nsims, xi$est, xi$vcov)
    l[[i]] <- list(sample = samp, loc.index = xi$loc.index, anc1.index = xi$anc1.index,
                   anc2.index = xi$anc2.index)
  }
  return(l)
}


#' Sample stratified survival parameters 
#'
#' Sample survival paramters stratified by variable. 
#' 
#' @param nsims Size of the posterior sample.
#' @param x A list of survival parameters. For example, see 'Treatment duration' in \link{sample_pars}.
#' 
#' @return A list of lists. The top level list is the value of the variable used for 
#' stratification. The second level list consists of matrices each containing random draws of the
#'  parameters of the survival distribution. The matrices are for the location parameter and
#'  each of the ancillary parameters. 
#'  
#' @keywords internal
sample_stratified_survpars <- function(nsims, x){
  l <- vector(length(x), mode = "list")
  names(l) <- names(x)
  for (i in 1:length(x)){
    l[[i]] <- sample_survpars(nsims, x[[i]])
  }
  return(l)
}

#' Sample parameters for the LCGM
#'
#' Sample parameters for the latent class growth model used to simulate the progression
#' of HAQ over time developed by Norton (2014).
#' 
#' @return List containing posterior sample of paramaters from the Norton(2014)
#' LCGM.
#' 
#' @export
sample_pars_haq_lcgm <- function(nsims){
  samp <- MASS::mvrnorm(nsims, haq.lcgm.pars$coef$coef, haq.lcgm.pars$vcov)
  if (class(samp) == "numeric") samp <- t(as.matrix(samp)) 
  lsamp <- list()
  indx.beta <- c()
  indx.delta <- unlist(haq.lcgm.pars$index[c("delta2", "delta3", "delta4")])
  indx.beta <- unlist(haq.lcgm.pars$index[paste0("beta", seq(1, 4))])
  lsamp[["delta"]] <- aperm(array(c(t(samp[, indx.delta])),
                                 dim = c(8, 3, nsims)),
                           perm = c(2, 1, 3))
  lsamp[["beta"]] <- aperm(array(c(t(samp[, indx.beta])),
                                  dim = c(4, 4, nsims)),
                            perm = c(2, 1, 3))
  return(lsamp)
}

#' Sample parameters for mixture model utility mapping
#'
#' Sample parameters for mixture model utility mapping. 
#' 
#' @return List containing posterior sample of paramaters from Hernandez Alva (2013) 
#' mixture model
#' 
#' @export
sample_pars_utility_mixture <- function(nsims, pain){
  samp <- MASS::mvrnorm(nsims, iviRA::util.mixture.pars$coef$coef, iviRA::util.mixture.pars$vcov)
  if (class(samp) == "numeric") samp <- t(as.matrix(samp)) 
  lsamp <- list()
  names <- names(util.mixture.pars$index)
  names <- names[!names %in% "delta"]
  for (n in names){
    indx <- util.mixture.pars$index[[n]]
    if (n %in% c("beta1", "beta2", "beta3", "beta4")){
      lsamp[[n]] <- samp[, indx, drop = FALSE] 
    } else{
      lsamp[[n]] <- samp[, indx]
    }
  }
  var.pars <- c(paste0("epsilon", seq(1, 4)), "mu")
  for (v in var.pars){
    lsamp[[v]] <- ifelse(lsamp[[v]] <= 0, 0, lsamp[[v]])
    lsamp[[v]] <- sqrt(lsamp[[v]])
  }
  indx.delta <- util.mixture.pars$index[["delta"]]
  lsamp[["delta"]] <- aperm(array(c(t(samp[, indx.delta])),
                                  dim = c(4, 3, nsims)),
                            perm = c(2, 1, 3))
  lsamp[["pain"]] <- pain
  return(lsamp)
}

#' Method of moments for Gamma Distribution
#'
#' Method of moments for Gamma Distribution
#' 
#' @return parameters of gamma distribution. Optionally return quantiles of
#' distribution as well.
#' 
#' @export
mom_gamma <- function(mean, sd, quant = NULL){
  scale <- sd^2/mean
  shape <- mean/scale
  l <- list(scale = scale, shape = shape)
  if(!is.null(quant)){
    q <- vector(length(quant), mode = "list")
    for (i in 1:length(q)){
      q[[i]] <- qgamma(quant[i], shape = shape, scale = scale)
    }
      names(q) <-  paste0("q", quant)
      l <- c(l, q)
  }
  return(l)
}

#' Method of moments for Beta Distribution
#'
#' Estimate parameters of Beta distribution from mean
#' and variance.
#' 
#' @return parameters of beta distribution
#' 
#' @export
mom_beta <- function(mean, sd, quant = NULL){
  alpha <- ((1 - mean) / sd^2 - 1 / mean) * mean ^ 2
  beta <- alpha * (1 / mean - 1)
  if(sum(alpha < 0) | sum(beta < 0)){
    stop(paste0("There is no Beta distribution with that mean and variance for
      at least one element in vector: returns
         alpha = [", paste(alpha, collapse = ", "), 
                "] and beta = [", paste(beta, collapse = ", "), "]"))
  }
  return(list(alpha = alpha, beta = beta))
}

#' Sample from gamma distributions
#'
#' Sample from multiple independent gamma distributions. The parameters of the gamma distribution are derived using methods of moments from
#' the means and standard errors. Produces a unique sample of size \code{n} using each element in 
#' \code{mean} and \code{sd}. 
#' 
#' @param n Number of samples to draw.
#' @param mean Vector of means.
#' @param se Vector of standard errors. 
#' @parm col_names Column names to returned matrix.
#' 
#' @return Matrix with each column a sample from a gamma distribution. 
#' 
#' @export
sample_gammas <- function(n, mean, se, col_names = NULL){
  which.fixed <- which(se == 0); which.gamma <- which(se > 0)
  gamma.pars <- mom_gamma(mean[which.gamma], se[which.gamma])
  samp.gamma <- rgamma(n * length(which.gamma), 
                 shape = gamma.pars$shape, scale = gamma.pars$scale)
  samp.gamma <- matrix(samp.gamma, nrow = n, ncol = length(which.gamma), byrow = T)
  samp.fixed <- matrix(mean[which.fixed], nrow = n, ncol = length(which.fixed), byrow = T)
  samp <- matrix(NA, nrow = n, ncol = length(se))
  samp[, which.gamma] <- samp.gamma
  samp[, which.fixed] <- samp.fixed
  colnames(samp) <- col_names
  return(samp)
}

#' Calculate mean, SD, and quantiles
#'
#' Calculate mean, SD, and quantiels of matrix or vector
#' 
#' @return matrix
#' @keywords internal
#' @export
apply_summary <- function(x){
  if (is.matrix(x)){
    y <- cbind(apply(x, 2, mean, na.rm = TRUE), apply(x, 2, sd, na.rm = TRUE),
               apply(x, 2, quantile, .025, na.rm = TRUE),
               apply(x, 2, quantile, .975, na.rm = TRUE))
  } else if (is.vector(x)) {
    y <- c(mean(x), sd(x), quantile(x, .025, na.rm = TRUE), quantile(x, .975, na.rm = TRUE))
    y <- t(matrix(y))
  } else if (is.array(x)){
     y <- cbind(c(t(apply(x, c(1, 2), mean))),
                c(t(apply(x, c(1, 2), sd))),
                c(t(apply(x, c(1, 2), quantile, .025))),
                c(t(apply(x, c(1, 2), quantile, .975))))
  }
  colnames(y) <- c("Mean", "SD", "Lower", "Upper")
  return(y)
}

#' Survival model parameters
#'
#' Summarize parameters from a survival model
#' 
#' @keywords internal
#' @return matrix
#' 
#' @export
survival_summary <- function(x){
  dists <- names(x)
  l <- vector(length(dists), mode = "list")
  n <- rep(NA, length(dists))
  for (i in 1:length(dists)){
    l[[i]] <- apply_summary(x[[dists[i]]]$sample)
    n[i] <- nrow(l[[i]])
  }
  l <- do.call(rbind, l)
  par.names <- paste0(rep(dists, times = n), " - ", rownames(l))
  return(list(par.summry = l, par.names = par.names))
}

#' Summary of parameters
#'
#' Table summarizing distribution of parameters
#' 
#' @return data.table
#' 
#' @export
par_table <- function(x, pat){
  acr.cats <- c("ACR <20", "ACR 20-50", "ACR 50-70", "ACR 70+")
  
  # rebound factor
  rebound <- data.table(Group = "Rebound factor", 
                        Distribution = "Uniform",
                        Parameter = "Increase in HAQ as a fraction of initial response",
                        Source = "Expert opinion")
  rebound <- cbind(rebound, apply_summary(x$rebound))
  
  # acr2eular
  n <- dim(x$acr2eular)[3]
  acr2eular.mat <- aperm(x$acr2eular, c(3,2,1))
  dim(acr2eular.mat) <- c(n, nrow(x$acr2eular) * ncol(x$acr2eular))
  acr2eular.dt <- data.table(Group = "Probability of EULAR response given ACR response", 
                        Distribution = "Dirichlet",
                        Parameter = c("ACR <20 to EULAR no response",
                                      "ACR <20 to EULAR moderate response",
                                      "ACR <20 to EULAR good response",
                                      "ACR 20-50 to EULAR no response",
                                      "ACR 20-50 to EULAR moderate response",
                                      "ACR 20-50 to EULAR good response",
                                      "ACR 50-70 to EULAR no response",
                                      "ACR 50-70 to EULAR moderate response",
                                      "ACR 50-70 to EULAR good response",
                                      "ACR 70+ to EULAR no response",
                                      "ACR 70+ to EULAR moderate response",
                                      "ACR 70+ to EULAR good response"),
                        Source = "Stevenson2016")
  acr2eular.dt <- cbind(acr2eular.dt, apply_summary(acr2eular.mat))
  
  # ACR to HAQ
  acr2haq.dt <- data.table(Group = "HAQ change by ACR response",
                        Distribution = "Normal",
                        Parameter = acr.cats,
                        Source = "Carlson2015")
  acr2haq.dt <- cbind(acr2haq.dt, apply_summary(x$acr2haq)) 
  
  # ACR to DAS28
  acr2das28.dt <- data.table(Group = "Change in DAS28 by ACR response",
                          Distribution = "Uniform",
                          Parameter = acr.cats,
                          Source = "Aletaha2005")
  acr2das28.dt <- cbind(acr2das28.dt, apply_summary(x$acr2das28)) 
  
  # ACR to SDAI
  acr2sdai.dt <- data.table(Group = "Change in SDAI by ACR response",
                          Distribution = "Uniform",
                          Parameter = acr.cats,
                          Source = "Aletaha2005")
  acr2sdai.dt <- cbind(acr2sdai.dt, apply_summary(x$acr2sdai)) 
  
  # ACR to CDAI
  acr2cdai.dt <- data.table(Group = "Change in CDAI by ACR response",
                         Distribution = "Uniform",
                         Parameter = acr.cats,
                         Source = "Aletaha2005")
  acr2cdai.dt <- cbind(acr2cdai.dt, apply_summary(x$acr2cdai)) 
  
  # EULAR to HAQ
  eular2haq.dt <- data.table(Group = "HAQ change by EULAR response",
                          Distribution = "Normal",
                          Parameter = c("No response", "Moderate response", "Good response"),
                          Source = "Stevenson2016")
  eular2haq.dt <- cbind(eular2haq.dt, apply_summary(x$eular2haq))
  
  # NMA ACR parameters
  ## coefficients
  acr.pars <- apply_summary(x$acr$pars)
  colnames(acr.pars) <- c("Mean", "SD", "Lower", "Upper")
  acr <- data.table(Group = "NMA ACR response",
                    Distribution = "Multivariate normal",
                    Parameter = c("cDMARDs mean", "Cutpoint ACR20", 
                                  "Cutpoint ACR50", "Cutpoint ACR70",
                                  therapy.pars$info$mname),
                    Source = "NMA")
  acr <- cbind(acr, acr.pars)
  
  ## relative risk reduction
  acr.rr <- data.table(Group = "NMA ACR response", 
                        Distribution = "Uniform",
                        Parameter = "Reduction in ACR response probabilities",
                        Source = "Carlson2015")
  acr.rr <- cbind(acr.rr, apply_summary(x$acr$rr))
  
  # NMA DAS28 parameters
  ## coefficients
  das28.pars <- apply_summary(x$das28$pars)
  colnames(das28.pars) <- c("Mean", "SD", "Lower", "Upper")
  das28 <- data.table(Group = "NMA DAS28",
                    Distribution = "Multivariate normal",
                    Parameter = c("cDMARDs mean",
                                  therapy.pars$info$mname),
                    Source = "NMA")
  das28 <- cbind(das28, das28.pars)
  
  ## relative risk reduction
  das28.rr <- data.table(Group = "NMA DAS28", 
                       Distribution = "Uniform",
                       Parameter = "Reduction in DAS28 response",
                       Source = "Carlson2015")
  das28.rr <- cbind(das28.rr, apply_summary(x$das28$rr))
  
  # NMA haq parameters
  ## coefficients
  haq.pars <- apply_summary(x$haq$pars)
  colnames(haq.pars) <- c("Mean", "SD", "Lower", "Upper")
  haq <- data.table(Group = "NMA HAQ",
                      Distribution = "Multivariate normal",
                      Parameter = c("cDMARDs mean",
                                    therapy.pars$info$mname),
                      Source = "NMA")
  haq <- cbind(haq, haq.pars)
  
  ## relative risk reduction
  haq.rr <- data.table(Group = "NMA HAQ", 
                         Distribution = "Uniform",
                         Parameter = "Reduction in HAQ response",
                         Source = "Carlson2015")
  haq.rr <- cbind(haq.rr, apply_summary(x$haq$rr))
  
  # treatment duration overall
  ttd.all.summary <- survival_summary(x$ttd.all)
  ttd.all <- data.table(Group = "Treatment duration - all patients",
                        Distribution = "Multivariate normal",
                        Parameter = ttd.all.summary$par.names,
                        Source = "Strand2013")
  ttd.all <- cbind(ttd.all, ttd.all.summary$par.summry)
  
  # treatment duration by eular response
  # moderate response
  ttd.em.summary <- survival_summary(x$ttd.eular$moderate)
  ttd.em <- data.table(Group = "Treatment duration - moderate EULAR response",
                       Distribution = "Multivariate normal",
                       Parameter = ttd.em.summary$par.names,
                       Source = "Stevenson2016")
  ttd.em <- cbind(ttd.em, ttd.em.summary$par.summry)
  
  ## good response
  ttd.eg.summary <- survival_summary(x$ttd.eular$good)
  ttd.eg <- data.table(Group = "Treatment duration - good EULAR response",
                       Distribution = "Multivariate normal",
                       Parameter = ttd.eg.summary$par.names,
                       Source = "Stevenson2016")
  ttd.eg <- cbind(ttd.eg, ttd.eg.summary$par.summry)
  
  # treatment duration by disease activity
  ## remission
  ttd.da.summary <- survival_summary(x$ttd.da)
  ttd.da <- data.table(Group = "Treatment duration by disease activity",
                       Distribution = "Multivariate normal",
                       Parameter = ttd.da.summary$par.names,
                       Source = "Zhang2011")
  ttd.da <- cbind(ttd.da, ttd.da.summary$par.summry)
  
  # Constant linear haq progression by therapy
  lhpt <- data.table(Group = "Yearly HAQ progression by therapy",
                    Distribution = "Normal",
                    Parameter = therapy.pars$info$mname,
                    Source = "Wolfe2010")
  lhpt <- cbind(lhpt, apply_summary(x$haq.lprog.therapy))
  
  # Constant linear haq progression by age
  lhpa <- data.table(Group = "Yearly HAQ progression by age relative to overall",
                    Distribution = "Normal",
                    Parameter = rownames(haq.lprog.age),
                    Source = "Michaud2011")
  lhpa <- cbind(lhpa, apply_summary(x$haq.lprog.age))
  
  
  # HAQ progression latent class growth model
  ## delta
  haq.lcgm.delta <- apply_summary(x$haq.lcgm$delta)
  delta.names <- list()
  delta.vars <- c("Intercept", "Age", "Female", "DAS28", "Disease duration", "Rheumatoid factor", "ACR criteria", "IMDQ4")
  for (i in 2:4){
    delta.names[[i]] <- paste0("Probability of class ", i, " membership - ", delta.vars)
  }
  
  ## beta
  haq.lcgm.beta <- apply_summary(x$haq.lcgm$beta)
  beta.names <- list()
  beta.vars <- c("Intercept", "Linear", "Quadratic", "Cubic")
  for (i in 1:4){
      beta.names[[i]] <- paste0("Class ", i, " predictors - ", beta.vars)
  }
  
  ## labels
  haq.lcgm <- data.table(Group = "HAQ progression LCGM",
                     Distribution = "Normal",
                     Parameter = c(do.call("c", delta.names),
                                   do.call("c", beta.names)),
                     Source = "Norton2014")
  haq.lcgm <- cbind(haq.lcgm, rbind(haq.lcgm.delta, haq.lcgm.beta))
  
  # mortality - coefficient baseline HAQ (logit scale)
  lo.mort <- data.table(Group = "Log odds mortality", 
                        Distribution = "Normal",
                        Parameter = "Baseline HAQ",
                        Source = "Wolfe2003")
  lo.mort <- cbind(lo.mort, apply_summary(x$logor.mort))
  
  # mortality - log hazard ratio change in HAQ from baseline
  lhr.mort.months <- c("0 - 6", ">6 - 12", ">12 - 24", ">24 - 36", ">36")
  lhr.mort <- data.table(Group = "Log hazard ratio mortality", 
                         Distribution = "Normal",
                         Parameter = paste0("Change in HAQ from baseline months ", lhr.mort.months),
                         Source = "Michaud2012")
  lhr.mort <- cbind(lhr.mort, apply_summary(x$mort.loghr.haqdif))
  
  # lifetables
  lt.labs <- c(paste0("qx - female, age ", x$lt$female[, "age"]), paste0("qx - male, age ", x$lt$male[, "age"]))
  lt <- data.table(Group = "Lifetables", 
                   Distribution = "Fixed",
                   Parameter = lt.labs,
                   Source = "CDC2011")
  lt <- cbind(lt, apply_summary(t(rbind(x$lt$female[, "qx", drop = FALSE],
                                        x$lt$male[, "qx", drop = FALSE]))))
  lt$SD <- 0
  
  # treatment cost
  tc.pars <- rep(c("Cost first 6 months", "Annual cost 6+ months"), 
                   each = nrow(x$treat.cost))
  tc <- data.table(Group = "Treatment cost", 
                   Distribution = "Fixed",
                   Parameter = paste0(tc.pars, " - ", x$treat.cost$mname),
                   Source = "Label")
  
  ## non-weight based therapies
  tmp = copy(x$treat.cost)
  tmp[, init_treat_cost := init_infusion_cost + init_rx_cost]
  tmp[, ann_treat_cost := ann_infusion_cost + ann_rx_cost]
  
  ## weight-based therapy
  tmpw <- tmp[sname == "ifxmtx"]
  tmpw.mat <- as.matrix(tmpw[, .(init_wgt_slope, ann_wgt_slope, init_util,
                                 ann_util, strength, price)])
  init.dose <- pat[, "weight"] %*% t(tmpw.mat[, c("init_wgt_slope"), drop = FALSE])
  ann.dose <- pat[, "weight"] %*% t(tmpw.mat[, c("ann_wgt_slope"), drop = FALSE])
  init.wcost <- matrix(NA, nrow = nrow(init.dose), ncol = ncol(init.dose))
  ann.wcost <- matrix(NA, nrow = nrow(ann.dose), ncol = ncol(ann.dose))
  for (j in 1:ncol(init.dose)){
    init.wcost[, j] <-  tmpw.mat[j, "price"] * tmpw.mat[j, "init_util"] * 
      ceiling(init.dose[, j]/tmpw.mat[j, "strength"])
    ann.wcost[, j] <-  tmpw.mat[j, "price"] * tmpw.mat[j, "ann_util"] * 
      ceiling(ann.dose[, j]/tmpw.mat[j, "strength"])
  }
  init.wcost <- apply(init.wcost, 2, mean)
  ann.wcost <- apply(ann.wcost, 2, mean) 
  tmp[sname %in% c("ifxmtx"), init_treat_cost := init.wcost + init_treat_cost]
  tmp[sname %in% c("ifxmtx"), ann_treat_cost := ann.wcost + ann_treat_cost]
  
  ## making table
  tc2 <- data.table(Mean = c(as.matrix(tmp[, .(init_treat_cost, ann_treat_cost)])))
  tc2[, ':=' (SD = 0, Lower = Mean, Upper = Mean)]
  tc2[, SD := ifelse(is.na(Mean), NA, SD)]
  tc <- cbind(tc, tc2)
  
  # hospital cost
  haq.grps <- c("0 - <0.5", "0.5 - <1", "1 - <1.5", "1.5 - <2",
                "2 - <2.5", ">2.5")
  # number of hospital days
  hdays <- data.table(Group = "Days in hospital per year",
                      Distribution = "Gamma",
                      Parameter = paste0("HAQ: ", haq.grps),
                      Source = "Carlson2015")
  hdays <- cbind(hdays, apply_summary(x$hosp.cost$hosp.days))
    
  # cost per hospital day
  hcost <- data.table(Group = "Cost per day in hospital",
                      Distribution = "Gamma",
                      Parameter = paste0("HAQ: ", haq.grps),
                      Source = "Carlson2015")
  hcost <- cbind(hcost, apply_summary(x$hosp.cost$cost.pday))
  
  # general management cost
  mgmt.summary <- apply_summary(x$mgmt.cost)
  mgmt.dist <- ifelse(mgmt.summary[, "SD"] == 0, "Fixed", "Gamma")
  mgmt <- data.table(Group = "General management cost",
                     Distribution = mgmt.dist,
                     Parameter = c("Chest x-ray", "X-ray visit", 
                                   "Outpatient followup", "Mantoux tuberculin skin test"),
                     Source = "Claxton2016", mgmt.summary)
  
  # productivity loss
  prod.loss <- data.table(Group = "Productivity loss",
                      Distribution = "Normal",
                      Parameter = "Regression coefficient - HAQ",
                      Source = "Wolfe2005")
  prod.loss <- cbind(prod.loss, apply_summary(x$prod.loss))
  
  # serious infection
  ## rate
  ttsi.summary <- survival_summary(x$ttsi)
  ttsi <- data.table(Group = "Time to serious infection",
                       Distribution = "Multivariate normal",
                       Parameter = ttsi.summary$par.names,
                       Source = "Singh2011")
  ttsi <- cbind(ttsi, ttsi.summary$par.summry)
  
  ## cost per infection
  si.cost <- data.table(Group = "Serious infection cost",
                        Distribution = "Normal",
                        Parameter = c("Cost per pneumonia hospitalization"),
                        Source = "ICER2016")
  si.cost <- cbind(si.cost, apply_summary(x$si.cost)) 
  
  ## utility loss per infection
  si.ul <- data.table(Group = "Serious infection utility loss",
                        Distribution = "Uniform",
                        Parameter = c("One month loss in utility"),
                        Source = "Stevenson2016")
  si.ul <- cbind(si.ul, apply_summary(x$si.ul)) 
  
  # utility mixture model
  ## labels
  b.vars <- list()
  for (i in 1:4){
    b.vars[[i]] <- paste0("Class ", i, " predictors - ", c("Intercept", "HAQ", "HAQ^2", "Pain/100", "Age/10k", "Age/10k^2"))
  }
  d.vars <- list()
  for (i in 1:3){
    d.vars[[i]] <- paste0("Probability of class ", i, " membership - ", c("Intercept", c("Haq", "Pain/100", "Pain/100^2")))
  }
  util <- data.table(Group = "Utility mixture model",
                      Distribution = "Multivariate normal",
                      Parameter = c(do.call("c", b.vars), 
                        "Predictor - Male",
                        paste0("Variance class ", seq(1, 4)),
                        "Variance", do.call("c", d.vars)),
                     Source = "Hernandez2013")
  ## extract parameters
  ### predictors - classes
  b <- list()
  for (i in 1:4){
    b[[i]] <- apply_summary(cbind(x$mixture.utility[[paste0("alpha", i)]], 
                                             x$mixture.utility[[paste0("beta", i)]]))
  }
  b <- do.call("rbind", b)
  
  ### predictor - male
  alpha.male <- apply_summary(x$mixture.utility$alpha)
  
  ### variance
  v <- list()
  for (i in 1:4){
    v[[i]] <- apply_summary(x$mixture.utility[[paste0("epsilon", i)]])
  }
  v <- do.call("rbind", v)
  v <- rbind(v, apply_summary(x$mixture.utility$mu))
  
  ### probability - class membership
  d <- apply_summary(x$mixture.utility$delta)
  
  ## combine labels and parameters
  util.parsum <- rbind(b, alpha.male, v, d)
  util <- cbind(util, util.parsum)
  
  # utility model wailoo
  util.wailoo <- data.table(Group = "Utility model Wailoo",
                      Distribution = "Normal",
                      Parameter =  c("Intercept", "Age", "Disease duration", "Baseline HAQ",
                                      "Male", "Previous DMARDs", "Current HAQ"),
                      Source = "Wailoo2006")
  util.wailoo <- cbind(util.wailoo, apply_summary(x$wailoo.utility)) 
                    
  # table 
  table <- rbind(rebound, acr2eular.dt, acr2das28.dt, acr2sdai.dt, acr2cdai.dt,
                 acr2haq.dt, eular2haq.dt,
                 ttd.all, ttd.em, ttd.eg, ttd.da,
                 acr, acr.rr, das28, das28.rr, haq, haq.rr,
                 lhpt, lhpa, haq.lcgm, lo.mort, lhr.mort, lt,
                 tc, hdays, hcost, mgmt, prod.loss, ttsi, si.cost, si.ul, util, util.wailoo)
  table <- table[, .(Group, Parameter, Mean, SD, Lower, Upper, Distribution, Source)]
  setnames(table, colnames(table), c("Group", "Parameter", "Posterior mean", "Posterior SD", 
                                     "Posterior 2.5%", "Posterior 97.5%", "Distribution", "Source"))
  return(table)
}


