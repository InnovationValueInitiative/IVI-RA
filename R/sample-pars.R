#' Sample parameters
#'
#' Randomly draw parameters from their (joint) probability distribution.
#' 
#' @param n Size of the posterior sample.
#' @param input_data An object of class 'input_data' returned from \link{get_input_data}.
#' @param model_structures An object of class 'model_structures' returned from 
#' \link{select_model_structures}.
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
#' @param tx_cost Treatment cost matrix and treatment lookup in format of iviRA::tx.cost.
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
#' @param haq_lprog_tx_mean Point estimate of linear yearly HAQ progression rate by treatment.
#' @param haq_lprog_tx_se Standard error of linear yearly HAQ progression rate by treatment.
#' @param haq_lcgm_pars Parameters of LCGM for HAQ progression.
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
#' @param tx_attr_ug_lower Lower bound for utility gain from treatment attributes.
#' @param tx_attr_ug_upper Upper bound for utility gain from treatment attributes.
#' @param tx_attr_ug_names Names of treatment attributes to be returned in sampled matrix
#' \code{tx.attr.ug}.
#' @param tx_names Vector of treatment names.
#' @param incidence_female Incidence rates for females. A \code{data.table} that must contain a
#'  row for each single-year of age from 0 to 100, a variable \code{events} (the number of events
#'   (i.e. cases of rheumatoid arthritis) for each age), and the variable \code{person_years}(the 
#'   time at risk for an event). 
#' @param incidence_male Identical to \code{incidence_female} but for males.
#' @param util_mixture_pain Summary statistics for bivariate distribution of HAQ and pain. Format
#' should be the same as iviRA::pain. Currently, each element of the list must be of length 1.
#' 
#' @return List containing samples for the following model parameters:
#' 
#' \describe{
#'  \item{tx.cost}{\code{data.table} of treatment cost data.}
#'  \item{lt}{A list with two elements for consisting of two matrics, one for males and
#'  one for females. Each matrix contains three variables: \code{age}, \code{qx} 
#'  (probability of death) and \code{logit_qx} (the logit of the probability of death). 
#'  Importantly, there is a row for each single-year of age from 0 to 100, which is passed to
#'  the \link{sim_iviRA} function.}
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
#'   \item{haq.lprog.tx}{A matrix of sampled yearly linear change in HAQ by treatment.
#'    The matrix has one column for each treatment in \code{iviRA::treatments}.}
#'    \item{haq.lprog.age}{A matrix of sampled yearly linear change in HAQ by age. The matrix
#'    has three columns for age < 40, age 40-64, and age 65+.}
#'    \item{haq.lcgm}{A list of two elements containing parameters from the latent class growth
#'     model. The first element is \code{delta} which is a an array of sampled matrices with
#'      each matrix containing coefficients predicting class membership. Rows are classes and columns index
#'     variables. \code{beta} is similar to \code{delta}, but each matrix contains coefficients 
#'     predicting HAQ as a function of time using a quadratic polynomial model.
#'    }
#'    \item{utility.mixture}{A list containing samples of all parameters in the Hernandez Alva (2013) mixture model. See 'Sampled mixture model parameters'
#'    for details.}
#'    \item{utility.wailoo}{A matrix of sampled regression coefficients from the model mapping HAQ to EQ5D utility in Wailoo (2006). Variables are
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
#'    \item{tx.attr.ug}{Matrix of sampled values of utility gains. Each column is a different treatment
#'    attribute.}
#'    \item{incidence}{A list with two elements for consisting of two matrics, one for males and
#'  one for females. Each column is a random sample of the incidence rate for a given age. 
#'  Importantly, there is a column for each single-year of age from 0 to 100, which is passed to
#'  the \link{sim_iviRA} function.}
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
sample_pars <- function(n = 100, input_data,
                        rebound_lower = .7, rebound_upper = 1,
                       ltfemale = iviRA::lifetable.female, ltmale = iviRA::lifetable.male,
                       acr2eular_mat = iviRA::acr2eular,
                       tx_cost = iviRA::tx.cost,
                       mort_logor = iviRA::mort.or$logor, mort_logor_se = iviRA::mort.or$logor_se,
                       mort_loghr_haqdif = iviRA::mort.hr.haqdif$loghr,
                       mort_loghr_se_haqdif = iviRA::mort.hr.haqdif$loghr_se,
                       ttd_all = iviRA::ttd.all, ttd_da = iviRA::ttd.da, ttd_eular = iviRA::ttd.eular,
                       nma_acr_mean = iviRA::nma.acr.naive$mean,
                       nma_acr_vcov = iviRA::nma.acr.naive$vcov,
                       nma_acr_rr_lower = .75, nma_acr_rr_upper = .92,
                       nma_das28_mean = iviRA::nma.das28.naive$mean,
                       nma_das28_vcov = iviRA::nma.das28.naive$vcov,
                       nma_das28_rr_lower = .75, nma_das28_rr_upper = .92,
                       nma_haq_mean = iviRA::nma.haq.naive$mean,
                       nma_haq_vcov = iviRA::nma.haq.naive$vcov,
                       nma_haq_rr_lower = .75, nma_haq_rr_upper = .92,
                       haq_lprog_tx_mean = iviRA::haq.lprog$tx$est,
                       haq_lprog_tx_se = iviRA::haq.lprog$tx$se,
                       haq_lcgm_pars = iviRA::haq.lcgm,
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
                       haq_lprog_age_mean = iviRA::haq.lprog$diff.age$est,
                       haq_lprog_age_se = iviRA::haq.lprog$diff.age$se,
                       hosp_days_mean = iviRA::hosp.cost$days_mean,
                       hosp_days_se = iviRA::hosp.cost$days_se,
                       hosp_cost_mean = iviRA::hosp.cost$cost_pday_mean,
                       hosp_cost_se = iviRA::hosp.cost$cost_pday_se,
                       mgmt_cost_mean = iviRA::mgmt.cost$est,
                       mgmt_cost_se = iviRA::mgmt.cost$se,
                       pl_mean = iviRA::prod.loss$est, pl_se = iviRA::prod.loss$se,
                       ttsi = iviRA::ttsi,
                       si_cost = 5873, si_cost_range = .2,
                       si_ul = .156, si_ul_range = .2,
                       tx_attr_ug_lower = iviRA::tx.attr$utility.gain$lower,
                       tx_attr_ug_upper = iviRA::tx.attr$utility.gain$upper,
                       tx_attr_ug_names = iviRA::tx.attr$utility.gain$var,
                       tx_names = iviRA::treatments$sname,
                       incidence_female = iviRA::incidence.female,
                       incidence_male = iviRA::incidence.male,
                       util_mixture_pain = iviRA::pain){
  
  # input data and model structures
  if (!inherits(input_data, "input_data")){
    stop("The argument 'input_data' must be of class 'input_data'")
  }

  # number of covariates in NMA treatment by covariate interactions
  nma.acr.ncovs <- ncol(input_data$x.acr)
  nma.haq.ncovs <- ncol(input_data$x.haq)
  nma.das28.ncovs <- ncol(input_data$x.das28)
  
  # sampling
  acr.cats <- c("ACR <20", "ACR 20-50", "ACR 50-70", "ACR 70+")
  sim <- list()
  sim$n <- n
  sim$rebound <- runif(n, rebound_lower, rebound_upper)
  sim$lt <- lt_data(ltmale, ltfemale)
  sim$tx.cost <- c(tx_cost,
                      list(discount = sample_uniforms(n, lower = tx_cost$cost$discount_lower, 
                                      tx_cost$cost$discount_upper, 
                                      col_names = tx_cost$cost$sname)))
  sim$acr2eular <- sample_dirichlets(n, acr2eular_mat)
  sim$acr2sdai <- sample_uniforms(n, acr2sdai_lower, acr2sdai_upper, acr.cats)
  sim$acr2cdai <- sample_uniforms(n, acr2cdai_lower, acr2cdai_upper, acr.cats)
  sim$acr2das28 <- sample_uniforms(n, acr2das28_lower, acr2das28_upper, acr.cats)
  sim$logor.mort <- sample_mvnorm(n, mort_logor, mort_logor_se^2)
  sim$mort.loghr.haqdif <- sample_normals(n, mort_loghr_haqdif, mort_loghr_se_haqdif,
                                         col_names = paste0("month_", c("less6", "6to12", "12to24", "24to36", "36to48")))
  sim$ttd.all <- sample_survpars(n, ttd_all)
  sim$ttd.da <- sample_survpars(n, ttd_da)
  sim$ttd.eular <- sample_stratified_survpars(n, ttd_eular) # make depend on distribution
  sim$eular2haq <- sample_normals(n, eular2haq_mean, eular2haq_se,
                                    col_names = c("no_response", "moderate_response", 
                                                  "good_response"))
  sim$acr <- sample_nma_acr(n, nma_acr_mean, nma_acr_vcov, rr_lower = nma_acr_rr_lower,
                              rr_upper = nma_acr_rr_upper, ncovs = nma.acr.ncovs,
                              tx_names = tx_names) 
  sim$das28 <- sample_nma_lm(n, nma_das28_mean, nma_das28_vcov, rr_lower = nma_das28_rr_lower,
                               rr_upper = nma_das28_rr_upper, ncovs = nma.das28.ncovs,
                               tx_names = tx_names) 
  sim$haq <- sample_nma_lm(n, nma_haq_mean, nma_haq_vcov, rr_lower = nma_haq_rr_lower,
                             rr_upper = nma_haq_rr_upper, ncovs = nma.haq.ncovs,
                             tx_names = tx_names) 
  sim$acr2haq <- sample_normals(n, acr2haq_mean, acr2haq_se, acr.cats) 
  sim$haq.lprog.tx <- sample_normals(n, haq_lprog_tx_mean,
                                                  haq_lprog_tx_se)
  sim$haq.lprog.age <- sample_normals(n, haq_lprog_age_mean, haq_lprog_age_se,
                                     col_names =  c("age_less40", "age40to64", "age_65plus"))
  sim$haq.lcgm <- sample_pars_haq_lcgm(n, pars = haq_lcgm_pars) 
  sim$utility.mixture <- sample_pars_utility_mixture(n, util_mixture_pain)    
  sim$utility.wailoo <- sample_normals(n, iviRA::utility.wailoo$est, iviRA::utility.wailoo$se,
                                         col_names = utility.wailoo$var) 
  hosp.cost.names <- c("haq_less0.5", "haq0.5to1", "haq1to1.5", "haq1.5to2", "haq2to2.5", "haq2.5plus")
  sim$hosp.cost <- list(hosp.days = sample_gammas(n, mean =hosp_days_mean, se =hosp_days_se,
                                                 col_names = hosp.cost.names),
                     cost.pday = sample_gammas(n, mean = hosp_cost_mean, se = hosp_cost_se,
                                              col_names = hosp.cost.names))
  sim$mgmt.cost <- sample_gammas(n, mean = mgmt_cost_mean, se = mgmt_cost_se,
                                col_names = c("chest_xray", "xray_visit", "outpatient_followup", "tuberculin_test"))
  sim$prod.loss <- rnorm(n, pl_mean, pl_se)
  sim$ttsi <- sample_normals(n, mean = ttsi$lograte, sd = ttsi$lograte_se,
                             col_names = tx_names)
  sim$si.cost <- runif(n, si_cost * (1 - si_cost_range),  si_cost * (1 + si_cost_range))
  sim$si.ul <- runif(n, si_ul * (1 - si_ul_range), si_ul * (1 + si_ul_range))
  sim$tx.attr.ug <-  sample_uniforms(n, tx_attr_ug_lower, tx_attr_ug_upper,
                                     tx_attr_ug_names)
  sim$incidence <- list(female = sample_gammas(n, shape = incidence_female$events, 
                                 rate = incidence_female$person_years),
                        male = sample_gammas(n, shape = incidence_male$events, 
                                             rate = incidence_male$person_years))
  return(sim)
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
#' @return List containing posterior samples of changes in outcomes.
#' 
#' @export
sample_nma_lm <- function(nsims, m, vcov, rr_lower = 1, rr_upper = 1, ncovs = 1, tx_names = NULL){
  rr.sim <- runif(nsims, rr_lower, rr_upper)
  sim <- sample_mvnorm(nsims, m, vcov)
  d <- array(sim[, -c(1)], dim = c(nrow(sim), ncovs, ncol(sim[, -c(1)])))
  dimnames(d)[[3]] <- colnames(sim[, -c(1)])
  return(list(rr = rr.sim, A = sim[, 1], d = d))
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
#' @param ncovs Number of treatment by covariate interactions.
#' @return List containing posterior sample of ACR response for each therapy
#' 
#' @export
sample_nma_acr <- function(nsims, m, vcov, rr_lower = 1, rr_upper = 1, 
                           ncovs = 1, tx_names = NULL){
  rr.sim <- runif(nsims, rr_lower, rr_upper)
  sim <- sample_mvnorm(nsims, m, vcov)
  d <- array(sim[, -c(1:4)], dim = c(nrow(sim), ncovs, ncol(sim[, -c(1:4)])))
  dimnames(d)[[3]] <- colnames(sim[, -c(1:4)])
  return(list(rr = rr.sim, A = sim[, 1], z2 = sim[, 3], z3 = sim[, 4], d = d))
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
sample_pars_haq_lcgm <- function(nsims, pars){
  samp <- MASS::mvrnorm(nsims, pars$coef$est, pars$vcov)
  if (is.vector(samp)) samp <- t(as.matrix(samp)) 
  index.delta <- which(pars$coef$parameter %in% paste0("delta", seq(2, 4)))
  index.beta <- which(pars$coef$parameter %in% paste0("beta", seq(1, 4)))
  lsamp <- list()
  lsamp[["delta"]] <- aperm(array(c(t(samp[, index.delta])),
                                 dim = c(8, 3, nsims)),
                           perm = c(2, 1, 3))
  lsamp[["beta"]] <- aperm(array(c(t(samp[, index.beta])),
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
  # sample
  samp <- MASS::mvrnorm(nsims, iviRA::utility.mixture$coef$est,
                        iviRA::utility.mixture$vcov)
  if (is.vector(samp)) samp <- t(as.matrix(samp)) 
  
  # separate sampels into parameter sets
  ## parameters except delta
  lsamp <- list()
  names <- unique(iviRA::utility.mixture$coef$parameter)
  names <- names[!names %in% "delta"]
  for (n in names){
    indx <- which(iviRA::utility.mixture$coef$parameter %in% n)
    if (n %in% c("beta1", "beta2", "beta3", "beta4")){
      lsamp[[n]] <- samp[, indx, drop = FALSE] 
    } else{
      lsamp[[n]] <- samp[, indx]
    }
  }
  
  ### variance cannot be negative
  var.pars <- c(paste0("epsilon", seq(1, 4)), "mu")
  for (v in var.pars){
    lsamp[[v]] <- ifelse(lsamp[[v]] <= 0, 0, lsamp[[v]])
    lsamp[[v]] <- sqrt(lsamp[[v]])
  }

  ## delta
  delta.names <- paste0("delta", seq(1, 3))
  indx.delta <- which(iviRA::utility.mixture$coef$parameter %in% delta.names)
  lsamp[["delta"]] <- aperm(array(c(t(samp[, indx.delta])),
                                  dim = c(4, 3, nsims)),
                            perm = c(2, 1, 3))
  
  ## pain
  lsamp[["pain"]] <- pain
  
  # return
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
#' @param col_names Column names to returned matrix.
#' @return Matrix with each column a sample from a gamma distribution. 
#' 
#' @export
sample_gammas <- function(n, mean = NULL, se = NULL, shape = NULL, rate = NULL,
                          col_names = NULL){
  if ((!is.null(mean) | !is.null(se)) & (!is.null(shape) | !is.null(rate))){
    stop("Specify mean and and se or shape and rate but not both.")
  }
  if ((!is.null(mean) & is.null(se)) | !is.null(se) & is.null(mean)){
    stop("mean and se must both be specified.")
  }
  if ((!is.null(shape) & is.null(rate)) | !is.null(rate) & is.null(shape)){
    stop("shape and rate must both be specified.")
  }
  if (!is.null(mean)){
      which.fixed <- which(se == 0); which.gamma <- which(se > 0)
      gamma.pars <- mom_gamma(mean[which.gamma], se[which.gamma])
      samp.fixed <- matrix(mean[which.fixed], nrow = n, ncol = length(which.fixed), byrow = T)
      samp <- matrix(NA, nrow = n, ncol = length(se))
  } else{
      which.fixed <- which(shape == 0); which.gamma <- which(shape > 0)
      gamma.pars <- list(scale = 1/rate[which.gamma], shape = shape[which.gamma]) 
      samp.fixed <- matrix(shape[which.fixed], nrow = n, ncol = length(which.fixed), byrow = T)
      samp <- matrix(NA, nrow = n, ncol = length(rate))
  }
  samp.gamma <- rgamma(n * length(which.gamma), 
                 shape = gamma.pars$shape, scale = gamma.pars$scale)
  samp.gamma <- matrix(samp.gamma, nrow = n, ncol = length(which.gamma), byrow = T)
  samp[, which.gamma] <- samp.gamma
  samp[, which.fixed] <- samp.fixed
  colnames(samp) <- col_names
  return(samp)
}

