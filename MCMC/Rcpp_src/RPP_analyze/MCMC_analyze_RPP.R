library(Rcpp)
library(RcppArmadillo)
library(RcppGSL)
library(ggplot2)
library(PRP)
data("RPP_filtered")
attach(RPP_filtered)

Rcpp::sourceCpp("PRP_MCMC_with_PI.cpp")
Rcpp::sourceCpp("PRP_MCMC_pubbias.cpp")

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

## DC criterion default ts
rpp_pval_ref = sapply(1:nrow(RPP_filtered),function(x)
  prior_prp(beta=c(beta_orig[x], beta_rep[x]),
           se=c(se_orig[x], se_rep[x]))$pval)
rpp_pval_ref_PI<-sapply(1:nrow(RPP_filtered),function(x) 
  prior_prp(beta=c(beta_orig[x], beta_rep[x]),
            se=c(se_orig[x],  se_rep[x]),report_PI = T)$predictive_interval)

## DC publication bias ts
rpp_pval_pubbias_ref = sapply(1:nrow(RPP_filtered),function(x)
  prior_prp(beta=c(beta_orig[x], beta_rep[x]),
            se=c(se_orig[x], se_rep[x]),test="pub_bias")$pval) 

## Distinguishability default ts 
rpp_pval<-sapply(1:nrow(RPP_filtered),function(x)
  PRP_MCMC_with_PI(beta=c(beta_orig[x], beta_rep[x]),
            se=c(se_orig[x], se_rep[x]), k_vec=root_vec)$pval)

## Distinguishability publication bias ts
rpp_pval_pubbias<-sapply(1:nrow(RPP_filtered),function(x)
  PRP_MCMC(beta=c(beta_orig[x], beta_rep[x]),
              se=c(se_orig[x], se_rep[x]), k_vec=root_vec))

## default ts DC v.s. D
rpp_ts=data.frame(priorPRP=rep(c("original","MCMC"),each=nrow(RPP_filtered)),rep_pval=c(rpp_pval_ref,rpp_pval))
pdf(file = paste0("rpp_pval_default_ts.pdf"),width=6,height=5)
ggplot(rpp_ts,aes(x=rep_pval))+
  geom_histogram(aes(fill=priorPRP),color="white",bins=15,alpha=0.5,boundary=0,position = 'identity')+
  scale_x_continuous(limits = c(0, 1))+
  scale_fill_manual("Criterion",labels=c("Distinguishability","DC"),values=c("#69b3a2", "#404080"))+
  xlab(expression(p["prior"]))+
  theme_bw()
dev.off()


## publication bias ts DC v.s. D
rpp_1=data.frame(priorPRP=rep(c("original","MCMC"),each=nrow(RPP_filtered)),rep_pval=c(rpp_pval_pubbias_ref,rpp_pval_pubbias))

pdf(file = paste0("rpp_pval_pubbias_ts.pdf"),width=6,height=5)
ggplot(rpp_1,aes(x=rep_pval))+
  geom_histogram(aes(fill=priorPRP),color="white",bins=20,alpha=0.5,boundary=0,position = 'identity')+
  scale_x_continuous(limits = c(0, 1))+
  scale_fill_manual("Criterion",labels=c("DC","Distinguishability"),values=c("#69b3a2", "#404080"))+
  xlab(expression(p["prior"]))+
  theme_bw()
dev.off()


pdf(file = paste0("rpp_comp.pdf"),width=6,height=6.5)
plot(rpp_pval_ref,rpp_pval,xlab=expression(P[prior]~-~DC),ylab=expression(P[prior]~-~Distinguishability),col="blue")
abline(c(0,0),c(1,1))
dev.off()

## Distingsuihability default vs publication bias ts
rpp_d=data.frame(test_statistics=rep(c("default","MCMC"),each=nrow(RPP_filtered)),rep_pval=c(rpp_pval,rpp_pval_pubbias))

pdf(file = paste0("rpp_pval_default_vs_pubbias.pdf"),width=6,height=5)
ggplot(rpp_d,aes(x=rep_pval))+
  geom_histogram(aes(fill=test_statistics),color="white",bins=10,alpha=0.5,boundary=0,position = 'identity')+
  scale_x_continuous(limits = c(0, 1))+
  scale_fill_manual("test statistic",labels=c("default",expression(T[pb]==hat(beta)[rep]/hat(beta)[orig])),values=c("#69b3a2", "#404080"))+
  xlab(expression(p["prior"]~Distinguishability~Criterion))+
  theme_bw()
dev.off()


##### predictive interval plot
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


data_new = data.frame(cbind("Study" =RPP_filtered[,1],"orig" =RPP_filtered[,2] ,
                            "rep"=RPP_filtered[,4] ,"pval"=pval,"PIlow" =lo,"PIhigh" =hi))

data_new[,1]<-as.factor(data_new[,1])

data<- data_new[order(data_new$pval,decreasing = T),]

data$Study2<- factor(data$Study,as.character(data$Study))
data$significance = as.factor((data$pval>0.05)+0)

pdf(file = paste0("rpp_PI_MCMC.pdf"),width=7,height=8.5)
ggplot(data)+
  geom_pointrange(aes(x=Study2, y=orig, ymin=PIlow, ymax=PIhigh),color="darkblue",shape=3,fatten=2)+
  geom_hline(yintercept = 0, linetype=2)+
  geom_point(aes(x=Study2,y=rep,color=significance))+
  labs( x = "study ID", y = "estimated effect") +
  scale_color_manual(" ",labels = c(expression(p[prior]~-~distinguishability< 0.05), expression(p[prior]~-~distinguishability> 0.05)), values = c("red", "blue")) +
  coord_flip()+
  ylim(-0.5,3.5)+
  theme_bw()
dev.off()



