library(Rcpp)
library(RcppArmadillo)
library(RcppGSL)
Rcpp::sourceCpp("PRP_MCMC.cpp")


library(PRP)
data("RPP_filtered")
attach(RPP_filtered)

rpp_pval<-sapply(1:nrow(RPP_filtered),function(x)
  PRP_MCMC(beta=c(beta_orig[x], beta_rep[x]),
            se=c(se_orig[x], se_rep[x]), k_vec=c(0.273)))

data = cbind(RPP_filtered,"PRP_MCMC" = rpp_pval)
write.table(data,file="res.dat")