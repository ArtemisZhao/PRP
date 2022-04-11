library(Rcpp)
library(RcppArmadillo)
library(RcppGSL)
library(PRP)
Rcpp::sourceCpp("scripts/PRP_MCMC.cpp")

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

args = commandArgs(trailingOnly=TRUE)
filename = args[1]

bb = args[2]


d = read.table(filename)
attach(d)
N = dim(d)[1]

pvec = sapply(1:N, function(x) PRP_MCMC(beta=c(V2[x], V4[x]), se=c(V3[x], V5[x]), k=root_vec))

#pvec_ref = sapply(1:N,function(x) prior_prp(beta=c(V2[x], V4[x]), se=c(V3[x], V5[x]))$pval)

outfile = paste0("output/batch_2grp_bb_sd_",bb,".prp.out")
outd = cbind(1:N, pvec)
write(file=outfile, t(outd), ncol=2)

pdf(file = paste0("plots/hist_batch_2grp_",bb,".pdf"),width=6,height=5.5)
hist(pvec,xlab="PRP_MCMC",col="grey")
dev.off()

# pdf(file = paste0("plots/comp_batch_2grp_",bb,".pdf"),width=6,height=6.5)
# plot(pvec_ref,pvec,xlab="PRP_original",ylab="PRP_MCMC",col="blue",cex=0.1)
# abline(c(0,0),c(1,1))
# dev.off()

