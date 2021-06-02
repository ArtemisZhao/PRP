#' Prior Predictive Replication p-value Calculation
#'
#' Assessing the prior predictive distribution and calculating the replication
#' p-value based on it.
#'
#' @param beta A vector, containing the estimates in the original study and the
#' replication study.
#' @param se A vector, containing the standard errors of the estimates in the original
#' study and thereplication study.
#' @param r_vec A vector, defining the prior reproducible model. Each r value
#' corresponds to a probability.
#' @param test A string, determining which test statistics to use. If not sepecified,
#' the default two-sided one will be used.
#' @param report_PI A boolean, denoting whether the 95% predictive interval for
#' the estimates be reported or not. This is option is only valid for two-sided
#' test statistics.
#'
#' @return
#' A list with the following components:
#' \item{grid}{ The grid points for the hyperparameters.}
#' \item{test_statistics}{ The test statistics used in calculating the replication p-value.}
#' \item{pvalue}{ The resulting prior predictive replicaiton p-value.}
#' \item{predictive_interval}{The 95% predictive interval if required.}
#'
#' @export
#'
prior_prp<-function(beta,se,r_vec = c(0, 8e-4, 6e-3, 0.024),test="two_sided",report_PI=F){
  reslist<-list()
  beta_o<-beta[1]
  beta_r<-beta[2]

  se2_o<-se[1]^2

  se2_r<-se[2]^2

  ####new grid
  #chis1<-qchisq(c(0.25,0.5,0.75),df=1)
  #eta2_vec = beta_o^2/chis1
  eta2_vec = (se2_o+beta_o^2)/qchisq(c(0.25, 0.5, 0.75), df=1)
  #pv = c(1.0, 0.99, 0.975, 0.95)
  rv = r_vec

  make_grid <-function(eta2){
    grid = sapply(rv, function(x)  c(eta2*(1-x), eta2*x))
    return(t(grid))
  }

  grid = c()

  for (i in 1:length(eta2_vec)){
    grid = rbind(grid, make_grid(eta2_vec[i]))
  }
  reslist[["grid"]] = grid

  omg2 = grid[,1]
  phi2 = grid[,2]

  mean<-sapply(1:length(omg2),function(x)
    beta_o/(se2_o/omg2[x]+phi2[x]/omg2[x]+1))

  var<-sapply(1:length(omg2),function(x)
    (se2_o+phi2[x])*omg2[x]/(se2_o+phi2[x]+omg2[x])+phi2[x]+se2_r)


  pval<-sapply(1:length(mean),function(x) pnorm(beta_r, mean=mean[x],sd=sqrt(var[x])))

  wts = dnorm(beta_o, mean=0, sd=sqrt(omg2+se2_o+phi2))
  wts = wts/sum(wts)

  pval_wt<-wts%*%pval

  res = NA
  reslist["test_statistics"]=test
  if (test=="pub_bias"){
    mean_scale=mean/beta_o
    sd_scale=sqrt(var)/abs(beta_o)
    pval_new<-sapply(1:length(mean),function(x) pnorm(beta_r/beta_o, mean=mean_scale[x],sd=sd_scale[x]))
    res = wts%*%pval_new
    reslist[["pvalue"]]=res
  }
  if (test=="two_sided"){
    res = 2*min(pval_wt,1-pval_wt)
    reslist[["pvalue"]]=res
    if (report_PI){
      mixture_CDF_right<-Vectorize(function(t) wts%*%pnorm(t,mean=mean,sd=sqrt(var))-0.975)
      root_right<-uniroot(mixture_CDF_right,c(mean(mean),10))$root

      mixture_CDF_left<-Vectorize(function(t) wts%*%pnorm(t,mean=mean,sd=sqrt(var))-0.025)
      root_left<-uniroot(mixture_CDF_left,c(-10,mean(mean)))$root

      reslist[["predictive_interval"]]=c(PI_left_95p=root_left,PI_right_95p=root_right)
    }
  }
  return(reslist)
}
