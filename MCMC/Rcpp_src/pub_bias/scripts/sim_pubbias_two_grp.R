set.seed(123)


args = commandArgs(trailingOnly=TRUE)

p_thresh = 1.0

if(length(args)>=1){
	p_thresh = as.numeric(args[1])
}


beta = rnorm(1,0,0.5)

samplesize<-c(50, 50)


wipi<-function(p){
	if (p<p_thresh){return(1)}
	else{return(0)}
}

betafinal<-c()
sdfinal<-c()
rep=5000


for (k in 1:rep){
	betalist<-c()
	sdlist<-c()
	for (j in 1:2){

		cur_indnum<-samplesize[j]
    #count=0
		while (TRUE){

		  x = rnorm(cur_indnum, 0, 1) 
		  y = 1+ beta*x + rnorm(cur_indnum,0,1)
		  
		  est_p = summary(lm(y~x))$coefficient[2,c(1,2,4)]
		  
			pass = 1
			if(j==1)
				pass = wipi(est_p[3])
			if(pass == 1){  
				betalist<-c(betalist,est_p[1])
				sdlist<-c(sdlist,est_p[2])
				break
			}
			#count=count+1
		}
	}
	betafinal<-rbind(betafinal,betalist)
	sdfinal<-rbind(sdfinal,sdlist)
}



outd = cbind(rep(1:rep), betafinal[,1], sdfinal[,1], betafinal[,2], sdfinal[,2])

outfile = paste0("sim_data/pubbias_2grp_pthresh_",round(p_thresh,2),".dat")
write(file=outfile, t(outd), ncol=5)

