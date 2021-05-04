#' Posterior Predictive Replication p-value Calculation
#'
#' @param beta A vector, containing the estimates in the original study and the
#' replication study.
#' @param se A vector, containing the standard errors of the estimates in the original
#' study and thereplication study.
#' @param L A value, determining the times of repeating simulation.
#' @param r_vec A vector, defining the prior reproducible model. Each r value
#' corresponds to a probability.
#' @param test A function designed to calculate the test quantity.
#' @param print_test_dist A boolean, determining whether the simulated test statistics
#' value difference will be plot as histogram or not.
#'
#' @return
#' A list with the following components:
#' \item{grid}{ The grid points for the hyperparameters.}
#' \item{test_statistics}{ The test statistics used in calculating the replication p-value.}
#' \item{n_sim}{ The L value.}
#' \item{test_stats_dif}{ The difference between the simulated test statistics
#' quantity and the original value.}
#' \item{pvalue}{ The resulting posterior predictive replicaiton p-value.}
#'
#' @importFrom mvtnorm dmvnorm
#'
#' @export
#'
posterior_prp<-function(beta,se,L=1000,r_vec = c(0,8e-4, 6e-3, 0.024),test,print_test_dist=FALSE){
  res<-list()
  sd2=se^2
  m<-length(beta)  ###number of replicates

  #chis1<-qchisq(c(0.25,0.5,0.75),df=1)
  ## test1 normal
  #eta2_vec = c(mean(beta),min(beta),max(beta))^2

  ## test 2 weighted average
  wts = 1/sd2
  center = sum(wts*beta)/sum(wts)
  eta2_vec = c(center, center-sqrt(1/sum(wts)),center+sqrt(1/sum(wts)))^2

  ## test 3
  #eta2_vec = c(center)^2
  #res[["eta_grid"]] = eta2_vec
  rv = r_vec

  make_grid <-function(eta2){
    grid = sapply(rv, function(x)  c(eta2*(1-x), eta2*x))
    return(t(grid))
  }

  grid = c()

  for (i in 1:length(eta2_vec)){
    grid = rbind(grid, make_grid(eta2_vec[i]))
  }
  res[["grid"]]= grid
  omg2_list = grid[,1]
  phi2_list = grid[,2]

  wts<-c()
  count=0
  for (i in 1:length(omg2_list)){
    omg2<-omg2_list[i]
    phi2<-phi2_list[i]
    Sigma<-matrix(omg2,ncol=m,nrow=m)+diag(c(sd2+phi2),nrow=m)
    wtsi<-dmvnorm(beta, mean=rep(0,m), sigma=Sigma)
    wts<-c(wts,wtsi)
  }
  wts<-wts/sum(wts)


  dist_list<-c()
  dist_list2<-c()
  for (t in 1:L){
    k<-sample(1:length(omg2_list),1,prob=wts)

    phi2<-phi2_list[k]
    omg2<-omg2_list[k]

    barbeta_pos_var<-1/(1/omg2+sum(1/(sd2+phi2)))
    barbeta_pos_mean<-barbeta_pos_var*sum(beta/(sd2+phi2))
    barbeta<-rnorm(1,barbeta_pos_mean,sqrt(barbeta_pos_var))

    betanewjs<-c()

    for (j in 1:m){
      #print(c(i,j))
      if (phi2==0){
        betaj=barbeta
      }
      else{
        betaj_var<-1/(1/phi2+1/sd2[j])
        betaj_mean<-betaj_var*(barbeta/phi2+beta[j]/sd2[j])
        betaj<-rnorm(1,betaj_mean,sqrt(betaj_var))
        #print(betaj)
      }

      betanewj = betaj+rnorm(1,0,sqrt(sd2[j]))
      betanewjs<-c(betanewjs,betanewj)
    }

      test_sim<-test(betanewjs,sd2,phi2,m)
      test_orig<-test(beta,sd2,phi2,m)
      dist_list = c(dist_list, test_sim)
      dist_list2= c(dist_list2, (test_sim - test_orig))
      count = count + ( test_sim > test_orig )
  }
  if (print_test_dist){
    #print(length(dist_list))
    hist(dist_list2)
  }
  res[["test_statistics"]]=test
  res[["n_sim"]]=L
  ###
  res[["test_stats_dif"]]=dist_list2
  ###
  res[["pvalue"]]=count/L

  return( res)
}

