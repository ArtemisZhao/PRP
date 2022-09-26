

Rcpp::sourceCpp("~/Desktop/rep_new/PRP/MCMC/Rcpp_src/RPP_analyze/PRP_MCMC_with_PI.cpp")

load("~/Desktop/rep_new/RPP_filtered_with_n.RData")
library(PRP)
load(RPP_filtered)

### jeff leek PI
PI<-function(betao,no,nr,lower=TRUE){
  z975=1.96
  int = sqrt(1/(no-3)+1/(nr-3))
  if (lower){
    pi.lo = betao-z975*int
    return(pi.lo)
  }
  else{
    pi.hi = betao+z975*int
    return(pi.hi)
  }
}

lo=c()
hi=c()
for (x in 1:nrow(data_with_n)){
  lo = c(lo,PI(data_with_n$beta_orig[x],data_with_n$n.O[x],data_with_n$n.R[x]))
  hi = c(hi,PI(data_with_n$beta_orig[x],data_with_n$n.O[x],data_with_n$n.R[x],lower=F))
}

data_jeff = data.frame(cbind("Study" =data_with_n[,1],"orig" =data_with_n[,2] ,
                            "rep"=RPP_filtered[,4] ,"pval"=NA, "PIlow" =lo,"PIhigh" =hi))
data_jeff[,1]<-as.factor(data_jeff[,1])

data_jeff$significance = as.factor(((data_jeff$rep>data_jeff$PIlow) & (data_jeff$rep<data_jeff$PIhigh))+0)


## Distinguishability PI
pval=c()
lo=c()
hi=c()
for (x in 1:nrow(RPP_filtered)){
  res = PRP_MCMC_with_PI(beta=c(beta_orig[x], beta_rep[x]),
                         se=c(se_orig[x], se_rep[x]), k_vec=root_vec)
  
  pval = c(pval,res$pval)
  lo = c(lo,res$left_bound)
  hi = c(hi,res$right_bound)
}

data_D = data.frame(cbind("Study" =RPP_filtered[,1],"orig" =RPP_filtered[,2] ,
                            "rep"=RPP_filtered[,4] ,"pval"=pval,"PIlow" =lo,"PIhigh" =hi))

data_D$significance = as.factor((data_D$pval>0.05)+0)
data_D[,1]<-as.factor(data_D[,1])

data_D<- data_D[order(data_D$pval,decreasing = T),]

### DC PI
rpp_pval<-sapply(1:nrow(RPP_filtered),function(x) 
  prior_prp(beta=c(beta_orig[x], beta_rep[x]),se=c(se_orig[x],  se_rep[x]))$pvalue)
rpp_PI<-sapply(1:nrow(RPP_filtered),function(x) 
  prior_prp(beta=c(beta_orig[x], beta_rep[x]),se=c(se_orig[x],  se_rep[x]),report_PI = T)$predictive_interval)

rpp_PI<-t(rpp_PI)

data_dc<-data.frame(cbind("Study"=RPP_filtered[,1],"orig" =RPP_filtered[,2] ,
               "rep"=RPP_filtered[,4] ,"pval"=rpp_pval, "PIlow" =rpp_PI[,1],"PIhigh" =rpp_PI[,2]))
data_dc$significance = as.factor((data_dc$pval>0.05)+0)


#### 
data_new <-cbind(rbind(data_jeff,data_dc,data_D),type=rep(c("Patil et al PI","DC PI","Distinguishability PI"),each=nrow(data_dc)))

data_new$Study <-factor(data_new$Study, levels=data_D$Study)
data_new$type<-factor(data_new$type, levels=c("Distinguishability PI","Patil et al PI","DC PI"))

pdf(file = paste0("rpp_overlay.pdf"),width=9,height=8.5)
ggplot(data_new)+
  geom_pointrange(aes(x=Study, y=orig, ymin=PIlow, ymax=PIhigh),color="darkblue",shape=3,fatten=2)+
  geom_hline(yintercept = 0, linetype=2)+
  geom_point(aes(x=Study,y=rep,color=significance))+
  labs( x = "study ID", y = "estimated effect") +
  scale_color_manual(" ",labels = c("Outside PI", "Inside PI"), values = c("red", "blue")) +
  coord_flip()+
  facet_grid(cols=vars(type))+
  ylim(-0.6,3.5)+
  theme_bw()
dev.off()
