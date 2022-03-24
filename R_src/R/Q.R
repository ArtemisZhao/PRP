#' Q test statistics
#'
#' This function calculates the Q test quantity.
#'
#' @param beta The original or the simulated estimated effects.
#' @param se2 The squared standard errors of the estimated effects.
#' @param barbeta The estimated true underlying effect.
#' @param phi2 The value of the hyperparameter phi.
#' @param m The number of replications
#'
#' @return The Q test statistic value
#'
#' @export
#' @keywords internal
#'
Q <- function(beta,se2,barbeta,phi2,m){
  q_test = sum((beta - barbeta)^2 / (se2 + phi2))
  return(q_test)
}
