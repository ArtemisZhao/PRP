set.seed(123)


args = commandArgs(trailingOnly=TRUE)

p_thresh = 1.0

if(length(args)>=1){
	p_thresh = as.numeric(args[1])
}


or<-c(2/3)

samplesize<-c(150, 150)


wipi<-function(p){
	if (p<p_thresh){return(1)}
	else{return(0)}
}

betafinal<-c()
sdfinal<-c()
rep=5000


for (k in 1:rep){
	curor = or
	betalist<-c()
	sdlist<-c()
	for (j in 1:2){

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
			pass = 1
			if(j==1)
				pass = wipi(est_p[3])
			if(pass == 1){  
				betalist<-c(betalist,est_p[1])
				sdlist<-c(sdlist,est_p[2])
				break
			}
		}
	}
	betafinal<-rbind(betafinal,betalist)
	sdfinal<-rbind(sdfinal,sdlist)
}




outd = cbind(rep(1:rep), betafinal[,1], sdfinal[,1], betafinal[,2], sdfinal[,2])

outfile = paste0("sim_data/pubbias_2grp_pthresh_",round(p_thresh,2),".dat")
write(file=outfile, t(outd), ncol=5)

