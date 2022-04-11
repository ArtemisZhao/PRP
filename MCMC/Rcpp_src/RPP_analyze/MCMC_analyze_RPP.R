library(Rcpp)
library(RcppArmadillo)
library(RcppGSL)
Rcpp::sourceCpp("PRP_MCMC_with_PI.cpp")

t= seq(0.001,0.05,0.001)
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

library(PRP)
data("RPP_filtered")
attach(RPP_filtered)

rpp_pval_ref = sapply(1:nrow(RPP_filtered),function(x)
  prior_prp(beta=c(beta_orig[x], beta_rep[x]),
           se=c(se_orig[x], se_rep[x]))$pval)
  
rpp_pval<-sapply(1:nrow(RPP_filtered),function(x)
  PRP_MCMC_with_PI(beta=c(beta_orig[x], beta_rep[x]),
            se=c(se_orig[x], se_rep[x]), k_vec=root_vec)$pval)

data = cbind(RPP_filtered,"PRP_MCMC" = rpp_pval)

write.table(data,file="res.dat")

pdf(file = paste0("rpp_comp.pdf"),width=6,height=6.5)
plot(rpp_pval_ref,rpp_pval,xlab="PRP_original",ylab="PRP_MCMC",col="blue")
abline(c(0,0),c(1,1))
dev.off()