#' Select model structures
#'
#' Select the model structures to use in the IVI-RA individual patient simulation.
#' 
#' @param tx_ihaq Model structure relating treatment to HAQ during the first 6 months of 
#' treatment. Options, which are equivalent to H1-H3 in the documentation are:
#' \itemize{
#' \item{\code{acr-haq}}{ H1: Treatment -> ACR -> HAQ }
#' \item{\code{acr-eular-haq}}{ H2: Treatment -> ACR -> EULAR -> HAQ}
#' \item{\code{haq}}{ H3:Treatment -> HAQ}
#' }
#' @param tx_iswitch Model structure relating treatment to switching during the first
#' 6 months of treatment. Options, which are equivalent to S1-S6 in the documentation are:
#' \itemize{
#' \item{\code{acr-switch}}{ S1: Treatment -> ACR -> Switch}
#' \item{\code{acr-das28-switch}}{ S2: Treatment -> ACR -> DAS28 -> Switch}
#' \item{\code{acr-sdai-switch}}{ S3: Treatment -> ACR -> SDAI -> Switch}
#' \item{\code{acr-cdai-switch}}{ S4: Treatment -> ACR -> CDAI -> Switch}
#' \item{\code{das28-switch}}{ S5: Treatment -> DAS28 -> Switch}
#'  \item{\code{acr-eular-switch}}{ S6: Treatment -> ACR -> EULAR -> Switch}
#' 
#' }
#' @param cdmards_haq_model Model used for long-term HAQ progression. Options are:
#' \itemize{
#' \item{lcgm}{ Latent class growth model}
#' \item{linear}{ Constant linear HAQ progression}
#' }
#' If \code{lgcm} is chosen, then a latent class growth model is used for cDMARDs
#' and NBT but a constant annual rate is is assumed for all other therapies; otherwise 
#' a constant linear HAQ progression is assumed for all therapies including cDMARDs and NBT.
#' @param ttd_cause Cause of treatment discontinuation. Options are:
#' \itemize{
#' \item{all}{ Treatment discontinuation due to any cause.}
#' \item{si}{ Treatment discontinuation due to serious infections}
#' }
#' @param ttd_dist Distribution used to model time to treatment discontinuaton. Options are:
#' \itemize{
#' \item{exponential}{ Exponential}
#' \item{weibull}{ Weibull}
#' \item{gompertz}{ Gompertz}
#' \item{gamma}{ Gamma}
#' \item{llogis}{ Log-logistic}
#' \item{lnorm}{ Lognormal}
#' \item{gengamma}{ Generalized gamma}
#' }
#' @param utility_model Model used to estimate patient utility as a function of HAQ and patient
#' characteristics. Options are:
#' \itemize{
#' \item{mixture}{ Hernandez Alava (2013) mixutre model}
#' \item{wailoo}{ Wailoo 2006 logistc regression}
#' }
#' 
#' @examples 
#' mod.structs <- select_model_structures(tx_ihaq = c("acr-haq", "acr-eular-haq"),
#'                                       tx_iswitch = c("acr-switch", "acr-eular-switch"),
#'                                       cdmards_haq_model = c("lcgm", "linear"),
#'                                       ttd_cause = c("all", "si"),
#'                                       ttd_dist = c("gengamma", "exponential"),
#'                                       utility_model = c("mixture", "wailoo"))
#' print(mod.structs)                                       
#' 
#' @export
select_model_structures <- function(tx_ihaq = "acr-haq",
                                    tx_iswitch = "acr-switch",
                                    cdmards_haq_model = "lcgm", 
                                    ttd_cause = "all",
                                    ttd_dist = "exponential",
                                    utility_model = "mixture"){
  # 
  n <- vector(length = 6)
  n[1] <- length(tx_ihaq)
  n[2] <- length(tx_iswitch)
  n[3] <- length(cdmards_haq_model)
  n[4] <- length(ttd_cause)
  n[5] <- length(ttd_dist)
  n[6] <- length(utility_model)
  if (max(n) > 1){
    n.g1 <- n[n > 1]
    max.n.g1 <- max(n.g1)
    if (any(max.n.g1 - n.g1 > 0)){
      stop("Length of all vectors must be the same")
    }
    if (any(max.n.g1 - n > 0)){
      if (n[1] == 1){
        tx_ihaq <- rep(tx_ihaq, max.n.g1)
      }
      if (n[2] == 1){
        tx_iswitch <- rep(tx_iswitch, max.n.g1)
      }
      if (n[3] == 1){
        cdmards_haq_model <- rep(cdmards_haq_model, max.n.g1)
      }
      if (n[4] == 1){
        ttd_cause <- rep(ttd_cause, max.n.g1)
      }
      if (n[5] ==1){
        ttd_dist <- rep(ttd_dist, max.n.g1)
      } 
      if (n[6] ==1){
        utility_model <- rep(utility_model, max.n.g1)
      } 
    }
  }
  
  # are valid options selected?
  if (any(!tx_ihaq %in% c("acr-haq", "acr-eular-haq", "haq"))){
    stop("Values in 'tx_ihaq' must be 'acr-haq', 'acr-eular-haq' or 'haq'.")
  }
  
  if (any(!tx_iswitch %in% c("acr-switch", "acr-das28-switch",
                             "acr-sdai-switch", "acr-cdai-switch", 
                             "das28-switch", "acr-eular-switch"))){
    stop(paste0("Values in 'tx_iswitch' must be 'acr-switch', 'acr-das28-switch', 'acr-sdai-switch',",
                " 'cr-cdai-switch', 'das28-switch', or 'acr-eular-switch'."))
  }
  
  if (any(!cdmards_haq_model %in% c("lcgm", "linear"))){
    stop("Values in 'cdmards_haq_model' must be 'lcgm' or 'linear'.")
  } 
  
  if (any(!ttd_cause %in% c("all", "si"))){
    stop(paste0("Values in 'ttd_cause' must be 'all', or 'si'."))
  } 
  
  if (any(!ttd_dist %in% c("exponential", "weibull", "gompertz", 
                           "gamma", "llogis", "lnorm", "gengamma"))){
    stop(paste0("Values in 'ttd_dist' must be 'exponential', 'weibull', 'gompertz', 'gamma',", 
                " 'llogis', 'lnorm', or 'gengamma'."))
  } 
  
  if (any(!utility_model %in% c("mixture", "wailoo"))){
    stop("Values in 'utility_model' must be 'mixture' or 'wailoo'.")
  } 
  
  # are valid combinations of options selected?
  ## tx_ihaq = acr-haq
  val <- ifelse(tx_ihaq == "acr-haq" & tx_iswitch == "acr-eular-switch", 1, 0)
  if (any(val > 0)){
    stop("'tx_iswitch' option 'acr-eular-switch' cannot be used with 'tx_ihaq' option
         'acr-haq'.")
  }
  
  ## tx_ihaq = haq
  val <- ifelse(tx_ihaq == "haq" & tx_iswitch != "das28-switch", 1, 0)
  if (any(val == 1)){
    stop("When 'tx_ihaq' option 'haq' is selected, 'tx_iswitch' must equal 'das28-switch'.")
  }
  
  ## ttd_cause == si
  val <- ifelse(ttd_cause == "si" & ttd_dist != "exponential", 1, 0)
  if (any(val == 1)){
    stop("When 'ttd_cause' option 'si' is selected, 'ttd_dist' must equal 'exponential'.")
  }
  
  # return
  model.structure <- matrix(c(tx_ihaq, tx_iswitch, cdmards_haq_model, 
                              ttd_cause, ttd_dist, utility_model), ncol = 6)
  colnames(model.structure) <- c("tx_ihaq", "tx_iswitch", "cdmards_haq_model", 
                                 "ttd_cause", "ttd_dist", "utility_model")
  class(model.structure) <- "model_structures"
  return(model.structure)
  }