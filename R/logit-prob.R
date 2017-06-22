#' Bayesian predicted probabilities for ordered logistic regression
#'
#' Posterior distribution of predicted probabilitities from an ordered logistic regression 
#' 
#' @param x Matrix containing variables in the model (i.e. the design matrix)
#'
#' @param beta Matrix of coefficients where each row contains random draws of the coefficients
#' from their posterior distribution
#'
#' @param cut Matrix of cutpoints, one for each category, where each row contains random draws
#' from their posterior distribution
#'
#' @return Matrix of predicted probabilities for each individual and simulation.
#'
#' @export
ologit_prob <- function(x, beta, cut){
  return(ologit_probC(x, beta, cut))
}