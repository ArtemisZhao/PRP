library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
filename = args[1]
bb = args[2]
#bb = 2
d = read.table(filename)

d = data.frame(d)
names(d)<-c("V1", "rep_pval")

pdf(file = paste0("plots/batch_2grp_",bb,".pdf"),width=4.5,height=4)
ggplot(data=d,aes(x=rep_pval))+
  geom_histogram(color="white",fill="#6186ad",bins=10,boundary=0,alpha=0.8)+
  scale_x_continuous(limits = c(0, 1))+
  theme_bw()+
  theme(plot.title = element_text(size=12,hjust = 0.5))+
  xlab(expression(p["prior"]))+
  ylim(0,3000)
dev.off()


############
sensitivity = c()
for (bb in seq(0,2,0.5)){
  d<-read.table(file=paste0("output/batch_2grp_bb_sd_",bb,".prp.out"))
  d = data.frame(d)
  names(d)<-c("V1", "rep_pval")
  sensitivity<-rbind(sensitivity,c(bb,length(which(d$rep_pval<0.05))/nrow(d)))
}
data<-data.frame(sensitivity)
colnames(data)<-c("eta","sensitivity")

pdf(file="plots/batch_2grp_sensitivity.pdf",width=4.5,height=4)
ggplot(data=data,aes(x=eta,y=sensitivity))+
  geom_line(color="#6186ad",size=0.8)+
  geom_point(color="#6186ad")+
  theme_bw()+
  xlab("batch effect magnitude")+
  ylim(0,1)+
  scale_x_continuous(breaks = seq(0, 2, by = 0.5))
dev.off()

