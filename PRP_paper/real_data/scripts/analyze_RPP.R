library(dplyr)
library(ggplot2)
library(PRP)
data("RPP_filtered")

attach(RPP_filtered)
rpp_pval<-sapply(1:nrow(RPP_filtered),function(x) prior_prp(beta=c(beta_orig[x], beta_rep[x]),se=c(se_orig[x],  se_rep[x]))$pvalue)
rpp_pval2<-sapply(1:nrow(RPP_filtered),function(x) prior_prp(beta=c(beta_orig[x], beta_rep[x]),se=c(se_orig[x],  se_rep[x]),test="pub_bias")$pvalue)

rpp_all=data.frame(test_statistics=rep(c("default","one sided"),each=nrow(RPP_filtered)),rep_pval=c(rpp_pval,rpp_pval2))

pdf(file = paste0("output/rpp_pval.pdf"),width=6,height=5)
ggplot(rpp_all,aes(x=rep_pval))+
  geom_histogram(aes(fill=test_statistics),color="white",bins=10,alpha=0.5,boundary=0,position = 'identity')+
  scale_x_continuous(limits = c(0, 1))+
  scale_fill_manual("test statistic",labels=c("default",expression(T[pb]==hat(beta)[rep]/hat(beta)[orig])),values=c("#69b3a2", "#404080"))+
  xlab(expression(p["prior"]))+
  theme_bw()
dev.off()

rpp_pval<-sapply(1:nrow(RPP_filtered),function(x) prior_prp(beta=c(beta_orig[x], beta_rep[x]),se=c(se_orig[x],  se_rep[x]))$pvalue)
rpp_PI<-sapply(1:nrow(RPP_filtered),function(x) prior_prp(beta=c(beta_orig[x], beta_rep[x]),se=c(se_orig[x],  se_rep[x]),report_PI = T)$predictive_interval)

rpp_PI<-t(rpp_PI)
rpp_PI<-cbind(Study, rpp_PI,rpp_pval, beta_orig,beta_rep)
rpp_PI<-data.frame(rpp_PI)

colnames(rpp_PI)<-c("Study","PIlow","PIhigh","pval","orig","observed")
rpp_PI[,1]<-as.factor(rpp_PI[,1])

data_new<- rpp_PI[order(rpp_PI$pval,decreasing = T),]

data_new$Study2<- factor(data_new$Study,as.character(data_new$Study))

data_new$significance = as.factor((data_new$pval>0.05)+0)


pdf(file = paste0("./output/rpp_PI.pdf"),width=6,height=8.5)
ggplot(data_new)+
  geom_pointrange(aes(x=Study2, y=orig, ymin=PIlow, ymax=PIhigh),color="darkblue",shape=3,fatten=2)+
  geom_hline(yintercept = 0, linetype=2)+
  geom_point(aes(x=Study2,y=observed,color=significance))+
  labs( x = "study ID", y = "estimated effect") +
  scale_color_manual(" ",labels = c(expression(p[prior]< 0.05), expression(p[prior]>0.05)), values = c("red", "blue")) +
  coord_flip()+
  ylim(-0.5,3.5)+
  theme_bw()
dev.off()

