#' Egger test statistics
#'
#' @param beta The original or the simulated estimated effects.
#' @param se2 The squared standard errors of the estimated effects.
#' @param barbeta The estimated true underlying effect.
#' @param phi2 The value of the hyperparameter phi.
#' @param m The number of replications
#'
#' @return The egger test statistic value
#'
#' @export
#'
egger <- function(beta,se2,barbeta,phi2,m){
  y = beta/sqrt(se2+phi2)
  x = 1/sqrt(se2+phi2)

  Sxx = sum( (x-mean(x))*x)
  Sxy = sum( (x-mean(x))*y)
  Syy = sum( (y-mean(y))*y)

  b1 = Sxy/Sxx
  b0 = mean(y)- b1*mean(x)

  s2 = (Syy - b1^2*Sxx)/m
  vb0 = s2*(1/m+mean(x)^2/Sxx)

  egger_test = b0^2/vb0
  return(egger_test)
}
