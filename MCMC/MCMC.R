t= seq(0.001, 0.05,0.001)
root_vec<-c(0)
for (t_ind in t){
  f = function(k){
    g = function(x){
      (1/(1+exp(0.5*(1/k)^2+x/k)))*dnorm(x)}
    return(integrate(g,-Inf,Inf)$value - t_ind) 
  }
  
  root<-uniroot(f,c(1e-100,1e5),tol = .Machine$double.eps^0.2)$root
  root_vec<-c(root_vec,root)
}

## independent MH algorithm

prior_PRP_mcmc <- function (beta, se, k_vec,nreps = 50000, burnin_rate = 0.01){
  hbo=beta[1]
  hbr=beta[2]
  so=se[1]
  sr=se[2]
  
  bbar = hbo
  k = sample(k_vec,1)
  
  bbar_sample = c(bbar)
  k_sample = c(k)
  k_max = k_vec[length(k_vec)]
  
  MH_step<-function( bbar, k, hbo, so ){
    
    # independent proposal
    k_p = sample(k_vec,1)
    bbar_p = rnorm(1, mean=hbo, sd = sqrt(k_max^2*hbo^2+so^2))
    
    
    #compute M-H ratio
    pr = dnorm(bbar_p, mean=hbo, sd = sqrt(k_max^2*hbo^2+so^2))/dnorm(bbar,  mean=hbo, sd = sqrt(k_max^2*hbo^2+so^2))
    lr = dnorm(hbo, mean = bbar_p, sd = sqrt(k_p^2*bbar_p^2+so^2))/dnorm(hbo, mean=bbar, sd=sqrt(k^2*bbar^2+so^2))
    mhr = lr/pr
    mhr = min(1, mhr)
    
    if(rbinom(1,1,mhr)){
      bbar = bbar_p
      k = k_p
    }
    
    return(c(bbar,k))
  }
  
  
  for (i in 1:nreps){
    rst = MH_step(bbar,k,hbo,so)
    bbar_sample = c(bbar_sample,rst[1])
    k_sample = c(k_sample, rst[2])
  }
  
  bbar_sample = bbar_sample[ceiling(nreps*burnin_rate):length(bbar_sample)]
  p1 = sapply(1:length(bbar_sample), function(x) pnorm(hbr, mean=bbar_sample[x], sd = sqrt(bbar_sample[x]^2*k_sample[x]^2+sr^2)))
  p2 = 1 - p1
  pval = 2*min(apply(cbind(p1,p2),2, mean))
  cat("p-value = ", pval,"\n")
  
  return(pval)
  
  #return(list(pval,bbar_vec))

}



library(PRP)
data("RPP_filtered")
attach(RPP_filtered)

res_new<-sapply(1,function(x) 
  prior_PRP_mcmc(beta=c(beta_orig[x], beta_rep[x]),
                 se=c(se_orig[x], se_rep[x]),k=root_vec, nreps = 50000))

# res<-c()
# for (t in 1:100){
# res_one=prior_PRP_mcmc(beta=c(beta_orig[1], beta_rep[1]),
#        se=c(se_orig[1], se_rep[1]),k=root_vec, nreps = 100000)
# res<-c(res,res_one)
# }
# mean(res)
# min(res)
# max(res)
# sd(res)


dat<-read.table("./res2.dat",header=TRUE)

newdata<-cbind(dat,"p_MCMC_new"=res_new)
#write.table(newdata,file="./res.dat")
#plot(newdata[,6],res_new,xlab="PRP_original",ylab="p_MCMC_new",col="blue")
#abline(c(0,0),c(1,1))

