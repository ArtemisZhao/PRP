set.seed(123)
args = commandArgs(trailingOnly=TRUE)


# key simulation parameters

rnse = 1.0  # replication noise level/se
bt_eff = 0.5 # beta true effect
bb_sd = 1.0 # beta batch effect, sd

bb_sd = as.numeric(args[1])

shuffle <- function(xv, rep){

	yv = xv
	for (i in 1:rep){
		index = sample(1:length(xv),2, replace=F)
		temp = yv[index[1]] 
		yv[index[1]] = yv[index[2]]
		yv[index[2]] = temp
	}
	return(yv)
}	



# simulate pure batch contamination
sim_batch<-function(bt=0, bb){

	y1 = bt*gv1+ bb*batch + rnorm(N)
	y2 = bt*gv2 + rnorm(N, sd = rnse)
	rst1 = summary(lm(y1~gv1))
	rst2 = summary(lm(y2~gv2))


	return(c(rst1$coef[2,1], rst1$coef[2,2], rst2$coef[2,1], rst2$coef[2,2]))
}



N = 100
gv1 = rbinom(N, 1, 0.4)
gv2 = rbinom(N, 1, 0.4)
batch = shuffle(gv1, rep=floor(N*0.25))


rst = t(sapply(1:5000, function(x) sim_batch(bt=bt_eff, bb=rnorm(1,sd=bb_sd))))



outd = cbind(1:5000, rst)
out_name = paste0("sim_data/batch_2grp_batch_sd_",round(bb_sd,2),".dat")
write(file= out_name, ncol=5, t(outd))



