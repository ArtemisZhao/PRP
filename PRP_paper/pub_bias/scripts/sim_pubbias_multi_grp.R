set.seed(123)





c = 0 # default, no censoring
args = commandArgs(trailingOnly=TRUE)
c = as.numeric(args[1]) 

####10 studies
studynum<-c(5,3,2)
####sample size in each study
indnum<-c(100,250,500) ###50 400 1000

samplesize = unlist(sapply(1:length(studynum), function(x) rep(indnum[x], studynum[x])))


bbar = log(2/3)
beta = bbar + rnorm(length(samplesize) ,sd = sqrt(0.02*bbar^2))
or = exp(beta)



#### censoring function
wipi<-function(p){
	return(exp(-c*p^(1.5)))
}




###without publication bias
# wipi<-function(p){
#     return(1)
#  }

betafinal<-c()
sdfinal<-c()
###repeats for 100 times
rep=5000
for (k in 1:rep){

	betalist<-c()
	sdlist<-c()
	for (j in 1:length(samplesize)){
		curor<-or[j]

		cur_indnum<-samplesize[j]
		beta_cur_est<-c()
		sd_cur_est<-c()

		while (TRUE){
			pcontrol<-runif(1,min=0.3,max=0.5)
			oddscontrol = pcontrol/(1-pcontrol) 
			oddstreat= oddscontrol * curor
			ptreat=oddstreat/(1+oddstreat)


			control<-rbinom(cur_indnum,1,pcontrol)
			treat<-rbinom(cur_indnum,1,ptreat)

			datay<-c(control,treat)
			datax<-c(rep(0,cur_indnum),rep(1,cur_indnum))
			est_p<-summary(glm(datay~datax,
					   family="binomial"))$coefficient[2,c(1,2,4)]
			if (rbinom(1,1,wipi(est_p[3])) == 1){
				betalist<-c(betalist,est_p[1])
				sdlist<-c(sdlist,est_p[2])
				break
			}
		}
	}


	betafinal<-rbind(betafinal,betalist)
	sdfinal<-rbind(sdfinal,sdlist)
}
outd = t(sapply(1:rep, function(x) c(x, as.vector(rbind(betafinal[x,], sdfinal[x,])))))
filename = paste0("sim_data/pubbias_multigrp_cparam_", c,".dat") 
write(file = filename, t(outd), ncol=dim(outd)[2])

