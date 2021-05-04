#' Q test statistics
#'
#' @param beta The original or the simulated estimated effects.
#' @param se2 The squared standard errors of the estimated effects.
#' @param phi2 The value of the hyperparameter phi.
#' @param m The number of replications
#'
#' @return The Q test statistic value
#' @export
#'
Q<-function(beta,se2,phi2,m){
  q_test = sum((beta - mean(beta))^2 / (se2 + phi2))
  return(q_test)
}
