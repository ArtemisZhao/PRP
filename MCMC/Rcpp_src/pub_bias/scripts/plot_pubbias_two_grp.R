library(ggplot2)

pubbias001<-read.table("./output/pubbias_2grp_pthresh_0.01.prp.out")
pubbias005<-read.table("./output/pubbias_2grp_pthresh_0.05.prp.out")
pubbiasno<-read.table("./output/pubbias_2grp_pthresh_1.prp.out")
pubbias001p<-read.table("./output/pubbias_2grp_pthresh_0.01_pubbias_ts.prp.out")
pubbias005p<-read.table("./output/pubbias_2grp_pthresh_0.05_pubbias_ts.prp.out")
pubbiasnop<-read.table("./output/pubbias_2grp_pthresh_1_pubbias_ts.prp.out")

data_plot<-data.frame(pval_censor=rep(c("0.01","0.05","no censor"),each=10000),
                      test_statistics=rep(c("default","pub_bias","default","pub_bias","default","pub_bias"),each=5000),
                      rep_pval=c(pubbias001$V2,pubbias001p$V2,pubbias005$V2,pubbias005p$V2,pubbiasno$V2,pubbiasnop$V2))

data_plot$pval_censor = factor(data_plot$pval_censor, levels=c('no censor','0.05','0.01'))
data_plot$test_statistics=factor(data_plot$test_statistics,levels=c("default","pub_bias"))
levels(data_plot$pval_censor)=c("'no censoring'", "'censoring at p=0.05'","'censoring at p=0.01'")
levels(data_plot$test_statistics) =c("'default test statistic'",expression(T[pb] == hat(beta)[rep]/hat(beta)[orig]))


pdf(file="./plots/mcmc_pub_bias.pdf", width=5.5, height=7)
ggplot(data=data_plot,aes(x=rep_pval))+
  geom_histogram(color="darkblue",fill="white",bins=10,boundary=1)+
  scale_x_continuous(limits=c(0,1))+
  facet_grid(pval_censor~test_statistics,labeller =label_parsed)+
  theme_bw()+xlab(expression(p[prior]~Distinguishability~Criterion))
dev.off()
