#' Treatment
#'
#' Dataset of treatments for rheumatoid arthritis.
#'
#' @format A data.table object with 1 row for each treatment and 4 columns:
#' \describe{
#'   \item{name}{Long-form name of treatment.}
#'   \item{mname}{Medium-form name of treatment.}
#'   \item{sname}{Short-form name of treatment.}
#'   \item{biologic}{Indicator equal to 1 if a treatment is a biologic and 0 otherwise.}
#' }
"treatments"

#' NMA estimates
#'
#' Estimates of ACR response, change in DAS28, and change in HAQ at 6 months for 
#' bDMARD naive patients from Bayesian random effects NMAs.
#'
#' @details ACR response parameters were estimated using an ordered probit model. Parameters for
#' DAS28 and HAQ were estimated using models for continuous data with a normal likelihood and 
#' identity link. The mean and variance-covariance matrix were estimated using the sample means
#' and sample covariance matrices from the posterior distributions of the Bayesian analyses. 
#' 
#' @format Each dataset is a list with two elements
#' \describe{
#'   \item{mean}{Mean of the parameters.}
#'   \item{vcov}{Variance-covariance matrix of the parameters.}
#' }
#' @name nma
"nma.acr.naive"
#' @rdname nma
"nma.das28.naive"
#' @rdname nma
"nma.haq.naive"

#' ACR to Eular Conversion
#'
#' Convert ACR response (ACR < 20, ACR20, ACR50, ACR70) to Eular response (none, moderate, good).
#'
#' @format A matrix where each row shows the probability that a particualar responder in an 
#' ACR category will have a given Eular response. Probabilities are based on observed 
#' relationships in the VARA database.

#' @source Stevenson, Matt, et al. "Adalimumab, etanercept, infliximab, certolizumab pegol, 
#' golimumab, tocilizumab and abatacept for the treatment of rheumatoid arthritis not previously
#'  treated with disease-modifying antirheumatic drugs and after the failure of conventional 
#'  disease-modifying antirheumatic drugs only: systematic review and economic evaluation." 
#' Health Technology Assessment 20.35 (2016): 1-610.
"acr2eular"

#' HAQ change by 6-month ACR Response
#'
#' Mean and standard deviation of HAQ change by 6-month ACR response category 
#' (ACR < 20, ACR20, ACR50, ACR70).
#'
#' @format A data.table with three columns. First column is ACR response category,
#' second column is mean change in HAQ at 6 months, and  third column is the standard error
#' of the mean change in HAQ.  
#' @source Carlson, Josh J., et al. "Economic evaluation of tocilizumab monotherapy compared to adalimumab monotherapy 
#' in the treatment of severe active rheumatoid arthritis." Value in Health 18.2 (2015): 173-179.
"acr2haq"

#' Change in SDAI by 6-month ACR Response
#'
#' Mean of change in SDAI by 6-month ACR response category (ACR < 20, ACR20, ACR50, ACR70).
#'
#' @format A list of two data.tables each with four columns. Column 1 is ACR response category,
#'  column 2 is the point estimate for the mean change in SDAI at 6 months, and columns 3/4 are 
#'  lower/upper bounds on the mean change in SDAI. Current default is to set bounds to +/-20% of the
#'   mean. The first list corresponds to evidence from the leflunomide dataset and the second list 
#'   to evidence from an inception cohort. 
#' @source Aletaha, Daniel, and Josef S. Smolen. "The simplified disease activity index (SDAI) and 
#' clinical disease activity index (CDAI) to monitor patients in standard clinical care." 
#' Best practice & research Clinical rheumatology 21.4 (2007): 663-675.
#' 
#' Smolen, J. S., et al. "A simplified disease activity index for rheumatoid arthritis for use in 
#' clinical practice." Rheumatology 42.2 (2003): 244-257.
#' 
#' Aletaha, Daniel, et al. "Acute phase reactants add little to composite disease activity indices
#' for rheumatoid arthritis: validation of a clinical activity score." Arthritis research & 
#' therapy 7.4 (2005): R796.
"acr2sdai"

#' Change in CDAI by 6-month ACR Response
#'
#' Mean of change in CDAI by 6-month ACR response category (ACR < 20, ACR20, ACR50, ACR70).
#'
#' @format A list of one data.tables with four columns. Column 1 is ACR response category,
#'  column 2 is the point estimate for the mean change in CDAI at 6 months, and columns 3/4 are 
#'  lower/upper bounds on the mean change in CDAI. Current default is to set bounds to +/-20% of
#'  the mean. Format is a list to facilitate inclusion of evidence from additional data sources 
#'  with different patient populations.
#' @source Aletaha, Daniel, and Josef S. Smolen. "The simplified disease activity index (SDAI) and 
#' clinical disease activity index (CDAI) to monitor patients in standard clinical care." 
#' Best practice & research Clinical rheumatology 21.4 (2007): 663-675.
#' 
#' Aletaha, Daniel, et al. "Acute phase reactants add little to composite disease activity indices
#' for rheumatoid arthritis: validation of a clinical activity score." Arthritis research & 
#' therapy 7.4 (2005): R796.
"acr2cdai"

#' Change in DAS28 by 6-month ACR Response
#'
#' Mean of change in DAS28 by 6-month ACR response category (ACR < 20, ACR20, ACR50, ACR70).
#'
#' @format A list of one data.table with four columns. Column 1 is ACR response category,
#'  column 2 is the point estimate for the mean change in DAS28 at 6 months, and columns 3/4 are 
#'  lower/upper bounds on the mean change in DAS28. Format is a list to facilitate inclusion of
#'   evidence from additional data sources with different patient populations.  
#' @source Aletaha, Daniel, and Josef S. Smolen. "The simplified disease activity index (SDAI) and 
#' clinical disease activity index (CDAI) to monitor patients in standard clinical care." 
#' Best practice & research Clinical rheumatology 21.4 (2007): 663-675.
#' 
#' Aletaha, Daniel, et al. "Acute phase reactants add little to composite disease activity indices
#' for rheumatoid arthritis: validation of a clinical activity score." Arthritis research & 
#' therapy 7.4 (2005): R796.
"acr2das28"

#' HAQ change by 6-month Eular Response
#'
#' Mean and standard deviation of HAQ change by 6-month Eular response category 
#' (none, moderate, good).
#'
#' @format A matrix with rows as Eular response categories. First column is mean and second 
#' column is standard deviation.
#' @source Stevenson, Matt, et al. "Adalimumab, etanercept, infliximab, certolizumab pegol, 
#' golimumab, tocilizumab and abatacept for the treatment of rheumatoid arthritis not previously
#'  treated with disease-modifying antirheumatic drugs and after the failure of conventional 
#'  disease-modifying antirheumatic drugs only: systematic review and economic evaluation." 
#' Health Technology Assessment 20.35 (2016): 1-610.
"eular2haq"


#' Linear HAQ progression rates
#'
#' Constant annual rate of HAQ progression by age and treatment.
#'
#' @format A list with two elements.
#' \describe{
#' \item{tx}{A matrix with 5 columns (sname, est, se, lower, upper) for the name of the treatment,
#' the annual progression rate for that treatment, the standard error of the progression,
#' the .025\% quantile of the progression rate and the 97.5\% quantile of the progression rate,
#' respectively.}
#' \item{diff.age}{A matrix with 5 columns (age, est, se, lower, upper) for the age band,
#' the difference between progression rate for a given age band and the overall rate,
#' the standard error of the progression rate difference, the .025\% quantile of the progression
#' rate difference, and the 97.5\% quantile of the progression rate difference, respectively. }
#' } 
#' 
#' @source Michaud, Kaleb, Gene Wallenstein, and Frederick Wolfe. "Treatment and nontreatment
#'  predictors of health assessment questionnaire disability progression in 
#'  rheumatoid arthritis: a longitudinal study of 18,485 patients." 
#'  Arthritis care & research 63.3 (2011): 366-372.
#'  
#' Wolfe, Frederick, and Kaleb Michaud. "The loss of health status in rheumatoid 
#' arthritis and the effect of biologic therapy: a longitudinal observational study." 
#' Arthritis research & therapy 12.2 (2010): 1.
"haq.lprog"

#' Latent Class Growth Model (LCGM) for HAQ progression
#'
#' Parameters for Norton (2014) LCGM used to simulate non-linear HAQ trajectories for
#' 4 latent classes.
#'
#' @format List containing coefficients (point estimates and parameter names)
#'  and variance-covariance matrix.
#' @source Norton, Sam, et al. "Health Assessment Questionnaire disability progression in early
#'  rheumatoid arthritis: systematic review and analysis of two inception cohorts." 
#'  Seminars in arthritis and rheumatism. Vol. 44. No. 2. WB Saunders, 2014.
"haq.lcgm"

#' Time to treatment discontinuation parameters
#'
#' \code{ttd.eular} and \code{ttd.da} are lists of lists containing treatment discontinuation
#' parameters stratified by EULAR response (moderate, high) and disease activity 
#' (remission, low, moderate, high) respectively; \code{ttd.all} is a list of lists containing
#' unstratified (i.e., all patients) treatment discontinuation parameters.
#' @details 
#' Time to treatment discontinuation paramters for each level of EULAR response or disease activity
#' are represented with a list as described in \link{sample_pars}. Models were fit using \emph{flexsurv}.
#' Survival curves are based on analyses of reconstructed individual patient data. 
#'
#' @format A list of lists of lists. The top level list is the level of EULAR response or
#'  disease activity. The second level list is the parametric distribution used. 
#'  The third level list contains the parameters of a given parametetric fit as described in
#'  \link{sample_pars}.

#' @source Stevenson, Matt, et al. "Adalimumab, etanercept, infliximab, certolizumab pegol, 
#' golimumab, tocilizumab and abatacept for the treatment of rheumatoid arthritis not previously
#'  treated with disease-modifying antirheumatic drugs and after the failure of conventional 
#'  disease-modifying antirheumatic drugs only: systematic review and economic evaluation." 
#' Health Technology Assessment 20.35 (2016): 1-610.
#' 
#' Zhang, Jie, et al. "Thresholds in disease activity for switching biologics in rheumatoid 
#' arthritis patients: experience from a large US cohort." 
#' Arthritis care & research 63.12 (2011): 1672-1679.
#' @name ttd
"ttd.eular"
#' @rdname ttd
"ttd.da"
#' @rdname ttd
"ttd.all"

#' Time to serious infection parameters
#'
#' Parameters for sampling serious infections using an exponential distribution.
#' @format A data.table with three columns:
#' \describe{
#' \item{sname}{Short name of treatment.}
#' \item{lograte}{Log of serious infection rate.}
#' \item{lograte_se}{Standard error of log of serious infection rate.}
#' }
#' @source Singh, Jasvinder A., et al. "Adverse effects of biologics: a network meta-analysis 
#' and Cochrane overview." Cochrane Database of Systematic Reviews 2 (2010).
"ttsi"

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

#' Incidence
#'
#' Annual incidence rate of rheumatoid arthritis in the United States.
#'
#' @name incidence
#' @format A data frame with 101 rows and 4 variables:
#' \describe{
#'   \item{age}{Age in years.}
#'   \item{events}{Number of events (i.e., number of patients with RA).}
#'   \item{person_years}{Person years, or, time at risk of event.}
#'   \item{incidence_rate}{Annual incidence rate.}
#'
#' }
#' @source Myasoedova, Elena, et al. "Is the incidence of rheumatoid arthritis rising?: 
#'results from Olmsted County, Minnesota, 1955–2007." 
#' Arthritis & Rheumatology 62.6 (2010): 1576-1582.
#' @name incidence
"incidence.female"
#' @rdname incidence
"incidence.male"

#' Mortality odds ratios by patient characteristics
#'
#' Impact of HAQ on odds ratio for mortality from table 4 in Wolfe et al (2003).
#'
#' @format A data frame with 1 row and 9 columns:
#' \describe{
#'   \item{var}{Name of variable}
#'   \item{or}{Odds ratio without radiographic data}
#'   \item{or_se}{Standard error of odds ratio without radiographic data}
#'   \item{or_lower}{Lower bound of 95 percent CI of odds ratio without radiographic data}
#'   \item{or_upper}{Upper bound of 95 percent CI of odds ratio without radiographic data}
#'   \item{logor}{Log odds ratio without radiographic data}
#'   \item{logor_lower}{Lower bound of 95 percent CI of log odds ratio without radiographic data}
#'   \item{logor_upper}{Upper bound of 95 percent CI of log odds ratio without radiographic data}
#'   \item{logor_se}{Standard error of log odds ratio without radiographic data}
#' }
#' 
#' @details The standard errors of the log odds ratios are derived from the 95 percent confidence 
#' interval for the log odds ratio, which is, in turn, calculated by taking the log of the lower and upper
#' bounds of the 95 perent confidence interval for the odds ratio. 
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

#' Treatment costs
#'
#' Drug acquisition and administration costs for medications used to treat rheumatoid
#' arthritis. 
#'
#' @format A data table with 5 variables.
#' \describe{
#'   \item{name}{Long-form name of treatment.}
#'   \item{sname}{Short-form name of treatment.}
#'   \item{dosage}{FDA approved dosage.}
#'   \item{strength_dosage_form}{Strength and dosage form of treatment}
#'   \item{init_dose_val}{Numeric value of dose approved by FDA during the first 6 months of
#'   treatment.}
#'   \item{ann_dose_val}{Numeric value of dose approved by FDA beyond the first 6 months.}
#'   \item{dose_unit}{Unit of FDA approved dose.}
#'   \item{init_num_doses}{Number of doses based on FDA used during the first 6 months.}
#'   \item{ann_num_doses}{Annual number of doses based on FDA beyond the first 6 months.}
#'   \item{strength_val}{Numeric value of strength of treatment.}
#'   \item{strength_unit}{Unit of treatment.}
#'   \item{wac_per_unit}{Wholesale acquisition cost per unit of treatment.}
#'   \item{infusion_cost}{Cost per infusion.}
#' }
"tx.cost"

#' Hospitalization costs
#'
#' Hospitalization costs due to RA.
#'
#' @format A data table with 5 variables.
#' \describe{
#'   \item{haq}{HAQ score range.}
#'   \item{days_mean}{Mean number of days in the hospital for a given HAQ score range.}
#'   \item{days_se}{Standard error of number of days in the hospital for a given HAQ score range.}
#'   \item{cost_pday_mean}{Mean number of costs per day in the hospital for a given HAQ score range.}
#'   \item{cost_pday_se}{Standard error of costs per day in the hospital for a given HAQ score range.}
#' }
#' @source Carlson, Josh J., et al. "Economic evaluation of tocilizumab monotherapy compared to adalimumab monotherapy 
#' in the treatment of severe active rheumatoid arthritis." Value in Health 18.2 (2015): 173-179.
"hosp.cost"

#' Cost of general management of RA
#'
#' Cost of general management of RA.
#'
#' @format A data table with 4 variables.
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

#' Productivity loss 
#'
#' Annual earnings loss from a 1-unit change in the HAQ score.
#'
#' @format A data table with 2 variables.
#' \describe{
#'   \item{est}{Estimate of the impact of 1-unit change in HAQ on annual lost earnings.}
#'   \item{se}{Standard error of point estimate of lost earnings.}
#' }
#' @source {Wolfe, Frederick, et al. "Household income and earnings losses among 6,396 persons
#'  with rheumatoid arthritis." The Journal of Rheumatology 32.10 (2005): 1875-1883.}
"prod.loss"

#' Mixture model utility mapping
#'
#' Coefficients and variance-covariance matrix from the Alava (2013) mixture model
#' used to map HAQ score to utilities.
#'
#' @format List containing coefficient vector (point estimates and parameter names) and 
#' variance-covariance matrix, and indices for parameters.
#'
#' @source Alava, Monica Hernandez, et al. "The relationship between EQ-5D, 
#'  HAQ and pain in patients with rheumatoid arthritis." Rheumatology 52.5 (2013): 944-950.
"utility.mixture"

#' Wailoo (2006) utility mapping
#'
#' Coefficients and standard errors from Wailoo (2006) logistic regression equation
#'  used to map HAQ score to utilities.
#'
#' @format A data.table containing the following columns.
#' \describe{
#'   \item{var}{Name of the variable.}
#'   \item{est}{Point estimate on the log odds scale.}
#'   \item{loghr_se}{Standard error of the point estimate.}
#' } 
#'
#' @source Wailoo, A., et al. "Modeling the cost effectiveness of etanercept, adalimumab and 
#' anakinra compared to infliximab in the treatment of patients with rheumatoid arthritis in the
#'  Medicare program." Rockville, MD: Agency for Healthcare Research and Quality (2006).
"utility.wailoo"

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

#' Treatment attributes
#'
#' Impact of treatment attributes other than safety and efficacy on utility.
#'
#' @format A list with two elements. 
#' \describe{
#'   \item{data}{A matrix where each column is a treatment attribute and each row corresponds
#'   to a treatment.}
#'   \item{utility.gain}{A matrix with two columns, lower and upper, corresponing to the lower
#'   and upper bounds of the utility gain associated with each treatment attribute in the 'data' 
#'   matrix.}
#' } 
#'
#' @source None. Currently determined by the user. 
"tx.attr"
