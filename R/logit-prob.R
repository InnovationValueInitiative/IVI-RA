#' Predicted probabilities for ordered logistic regression
#'
#' Predicted probabilitities from an ordered logistic regression for a single individual. This function
#' is currently not vectorized. 
#' 
#' @param x Row vector containing variables in the model for a single individual. 
#'
#' @param beta Vector of coefficients.
#'
#' @param cut Vector of cutpoints, one for each category. 
#'
#' @return Vector of predicted probabilities for the individual. 
#'
#' @export
ologit_prob <- function(x, beta, cut){
  return(ologit_probC(x, beta, cut))
}

#' Predicted probabilities for multinomial logistic regression
#'
#' Predicted probabilitities from a multinomial logistic regression for a single individual. This function
#' is currently not vectorized. 
#' 
#' @param x Row vector containing variables in the model for a single individual. 
#'
#' @param beta Matrix of coefficients; one for each outcome other than the reference outcome.
#'
#' @return Vector of predicted probabilities for the individual. 
#'
#' @export
mlogit_prob <- function(x, beta){
  return(mlogit_probC(x, beta))
}