#' US Lifetables 2011
#'
#' Life tables by single-year of age from National Vital Statistics Reports
#' Volume 64, Number 11.
#'
#' @name lifetable
#' @format A data frame with 101 rows and 7 variables:
#' \describe{
#'   \item{age}{Age in years.}
#'   \item{qx}{Probability of dying between ages x and x + 1.}
#'   \item{lx}{Number surviving to age x.}
#'   \item{dx}{Number dying between ages x and x + 1.}
#'   \item{L}{Person-years lived between ages x and x + 1.}
#'   \item{Tx}{Total number of person-years lived above age x.}
#'   \item{ex}{Expectation of life at age x.}
#'
#' }
#' @source \url{http://www.cdc.gov/nchs/data/nvsr/nvsr64/nvsr64_11.pdf}
#' @name lifetable
"lifetable.female"
#' @rdname lifetable
"lifetable.male"

#' Mortality odds ratios by patient characteristics
#'
#' Impact of HAQ on odds ratio for mortality from table 4 in Wolfe et al (2003).
#'
#' @format A data frame with 1 row and 9 columns:
#' \describe{
#'   \item{var}{Name of variable}
#'   \item{or}{Odds ratio without radiographic data}
#'   \item{or_se}{Standard error of odds ratio without radiographic data}
#'   \item{logor}{Log odds ratio without radiographic data}
#'   \item{logor_se}{Standard error of log odds ratio without radiographic data}
#' }
#' 
#' @details The standard errors of the log odds ratios are derived from the standard error of the 
#' odds ratios. By the delta method, \eqn{var(OR) = exp(\beta)var(\beta)exp(\beta)} where 
#' \eqn{\beta} is the regression coefficient of the logistic regression. Since \eqn{var(\beta)} 
#' is the variance of the log odds ratio, it follows that
#'  \eqn{var[log(OR)] = [var(OR)/exp(\beta)}]^2.
#' 
#' @source {Wolfe, Frederick, et al. "Predicting mortality in patients with rheumatoid arthritis." 
#' Arthritis & Rheumatism 48.6 (2003): 1530-1542.}
"mort.or"

#' Mortality hazard ratio for change in HAQ score
#'
#' Impact of change in HAQ on for mortality from table 3 in Michaud et al (2012).
#'
#' @format A data frame with 5 rows and 8 variables. Key variables are month, log of hazard ratio, 
#' and standard error of log of hazard ratio.
#' \describe{
#'   \item{month}{Month or reported hazard ratio}
#'   \item{loghr}{Log of hazard ratio for impact of .25 unit HAQ increase on mortality}
#'   \item{loghr_se}{Standard error of log of hazard ratio for impact of .25 unit HAQ increase on mortality}
#' }
#' @source {Michaud, Kaleb, Montserrat Vera-Llonch, and Gerry Oster. "Mortality risk by functional status and health-related quality of life in patients with rheumatoid arthritis." 
#' The Journal of rheumatology 39.1 (2012): 54-59.}
"mort.hr.haqdif"

#' Therapy parameters
#'
#' Parameters for effect of therapy on 6 month ACR response and haq progression rates. Also 
#' includes parameters needed to estimate treatment costs.
#'
#' @format A list of the following lists
#' \describe{
#'   \item{info}{Basic information about therapies in model}
#'   \item{icon_acr_nma}{Effect of treatment on ACR response probability from the icon NMA.
#'   Contains mean/sd for parameters in the NMA (Baseline effect for cDMARDs A, cutpoints z, and
#'  treatment effects delta). Also contains predicted probability of ACR category with parameters
#'  at posterior means.}
#'   \item{nice_acr_nma}{Effect of treatment on ACR response probability from the NICE NMA.
#'   Contains probability of each ACR response category by therapy as well as estimated number
#'   of observations needed to draw from a Dirichlet distribution. }
#'   \item{haq_prog}{Effect of treatment on HAQ progression rate. In particular, contains mean
#'    and standard error of yearly change in HAQ progression by therapy}
#'   \item{cost}{Treatment cost by therapy.}
#' }
#' @source {
#'   \describe{
#'   \item{nice_acr_nma}{Table 37, Stevenson, Matt, et al. "Adalimumab, etanercept, infliximab, certolizumab pegol, 
#'  golimumab, tocilizumab and abatacept for the treatment of rheumatoid arthritis not previously
#'  treated with disease-modifying antirheumatic drugs and after the failure of conventional 
#'  disease-modifying antirheumatic drugs only: systematic review and economic evaluation." 
#' Health Technology Assessment 20.35 (2016): 1-610.}
#' \item{haq_prog}{Wolfe, Frederick, and Kaleb Michaud. "The loss of health status in rheumatoid 
#' arthritis and the effect of biologic therapy: a longitudinal observational study." 
#' Arthritis research & therapy 12.2 (2010): 1.}
#' }}
"therapy.pars"

#' Survival analysis parameters for therapy duration by Eular response category
#'
#' A list of parameters estimated in a parametric survival analysis of treatment duration.
#' Models were fit using \emph{flexsurv}. Survival curves were scanned from the 2016 NICE  
#' HTA and then converted to individual-level patient dat using the \code{reconstruct_ipd}
#' function. Survival curves are stratified by patients with moderate and good responses 
#' to treatment according to Eular response categories.
#'
#' @format A list of lists. One list for each parametric distribution. Within each distribution
#' the list contains parameters produced by \code{flexsurv} on real-line scale
#'  (e.g. after log transformation if they are defined as positive). The parameters are:
#' \describe{
#'   \item{coef}{Coefficients for all parameters including location parameter and 
#'   ancillary parameter(s)}
#'   \item{vcov}{Covariance matrix for all parameters including location parameter and 
#'   ancillary parameter(s)}
#' }
#' @source Stevenson, Matt, et al. "Adalimumab, etanercept, infliximab, certolizumab pegol, 
#' golimumab, tocilizumab and abatacept for the treatment of rheumatoid arthritis not previously
#'  treated with disease-modifying antirheumatic drugs and after the failure of conventional 
#'  disease-modifying antirheumatic drugs only: systematic review and economic evaluation." 
#' Health Technology Assessment 20.35 (2016): 1-610.
#' @name dur.eular
"dur.eular.mod"
#' @rdname dur.eular
"dur.eular.good"
#' @rdname dur.eular
"dur.eular.mod.adj"
#' @rdname dur.eular
"dur.eular.good.adj"

#' ACR to Eular Conversion
#'
#' Convert ACR response (ACR < 20, ACR20, ACR50, ACR70) to 
#' Eular response (none, moderate, good)
#'
#' @format A list with two matrices. Each row in each matrix shows the probability that a 
#' particualar responder in an ACR category will have a given Eular response. Probabilities in
#' the first matrix ("observed") are based on observed relationships in the VARA database. In
#' the second matrix ("strict") the probabilities assume that patients with ACR < 20 have no
#' Eular response, that patients with ACR20 and ACR50 have moderate Eular response, and that 
#' patients with ACR70 have a good Eular response. 

#' @source Stevenson, Matt, et al. "Adalimumab, etanercept, infliximab, certolizumab pegol, 
#' golimumab, tocilizumab and abatacept for the treatment of rheumatoid arthritis not previously
#'  treated with disease-modifying antirheumatic drugs and after the failure of conventional 
#'  disease-modifying antirheumatic drugs only: systematic review and economic evaluation." 
#' Health Technology Assessment 20.35 (2016): 1-610.
"acr2eular"

#' HAQ change by 6-month Eular Response
#'
#' HAQ change by 6-month Eular response
#'
#' @format A matrix with rows as Eular response categories. First column is mean and second 
#' column is standard deviation.
#' @source Stevenson, Matt, et al. "Adalimumab, etanercept, infliximab, certolizumab pegol, 
#' golimumab, tocilizumab and abatacept for the treatment of rheumatoid arthritis not previously
#'  treated with disease-modifying antirheumatic drugs and after the failure of conventional 
#'  disease-modifying antirheumatic drugs only: systematic review and economic evaluation." 
#' Health Technology Assessment 20.35 (2016): 1-610.
"eular2haq"

#' Linear HAQ progression rate differences by age
#'
#' Impact of age on annual linear HAQ progression rate.
#'
#' @format A matrix with rows as age category. First column is mean and second 
#' column is standard error
#' @source Michaud, Kaleb, Gene Wallenstein, and Frederick Wolfe. "Treatment and nontreatment
#'  predictors of health assessment questionnaire disability progression in 
#'  rheumatoid arthritis: a longitudinal study of 18,485 patients." 
#'  Arthritis care & research 63.3 (2011): 366-372.
"haq.lprog.age"

#' Mixture model utility mapping
#'
#' Coefficients and variance-covariance matrix from the Alava (2013) mixture model
#' used to map HAQ score to utilities.
#'
#' @format List containing coefficient vector, variance-covariance matrix, and indices for
#' parameters.
#'
#' @source Alava, Mónica Hernández, et al. "The relationship between EQ-5D, 
#' HAQ and pain in patients with rheumatoid arthritis." Rheumatology 52.5 (2013): 944-950.
"util.mixture.pars"

#' Wailoo (2006) utility mapping
#'
#' Coefficients and standard errors from Wailoo (2006) model
#' used to map HAQ score to utilities.
#'
#' @format List containing vectors of coefficients and standard errors
#'
#' @source Wailoo, A., et al. "Modeling the cost effectiveness of etanercept, adalimumab and 
#' anakinra compared to infliximab in the treatment of patients with rheumatoid arthritis in the
#'  Medicare program." Rockville, MD: Agency for Healthcare Research and Quality (2006).
"util.wailoo.pars"

#' Correlation between HAQ and pain
#'
# Correlation between HAQ and pain as meausred on the visual analog scale (VAS).
#'
#' @format List containing mean and variance of HAQ and pain scores, as well as the 
#' correlation between pain and HAQ.
#'
#' @source Sarzi-Puttini, Piercarlo, et al. "Correlation of the score for subjective pain with
#'  physical disability, clinical and radiographic scores in recent onset rheumatoid arthritis."
#'  BMC musculoskeletal disorders 3.1 (2002): 1..
"pain"

#' Cost of general management of RA
#'
#' Cost of general management of RA
#'
#' @format A data frame with 4 variables.
#' \describe{
#'   \item{service}{Medical service performed}
#'   \item{est}{Point estimate for cost}
#'   \item{lower}{Lower range for cost}
#'   \item{upper}{Upper range for cost}
#' }
#' @source {Claxton, Lindsay, et al. "An Economic Evaluation of Tofacitinib Treatment in Rheumatoid 
#' Arthritis: Modeling the Cost of Treatment Strategies in the United States." Journal of managed
#'care & specialty pharmacy 22.9 (2016): 1088-1102.}
"mgmt.cost"


