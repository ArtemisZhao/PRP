library(PRP)



run_analysis<-function(x){

	p = length(x)
	beta = x[seq(2,p,2)]
	se = x[seq(3,p,2)]

	return(posterior_prp(beta,se,test=Q,r_vec=c(0))$pvalue)
}



d = read.table("sim_data/batch_multigrp_batch_sd_0.dat")
attach(d)
N = dim(d)[1]

pvec = apply(d, 1, function(x) run_analysis(x))

outfile = paste0("output/batch_multigrp_batch_sd_0_fixed_eff.prp.out")
outd = cbind(1:N, pvec)
write(file=outfile, t(outd), ncol=2)


