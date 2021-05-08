library(PRP)



run_analysis<-function(x){

	p = length(x)
	beta = x[seq(2,p,2)]
	se = x[seq(3,p,2)]

	return(posterior_prp(beta,se,test=Q)$pvalue)
}


args = commandArgs(trailingOnly=TRUE)
filename = args[1]
bb = args[2]


d = read.table(filename)
attach(d)
N = dim(d)[1]

pvec = apply(d, 1, function(x) run_analysis(x))

outfile = paste0("output/batch_multigrp_batch_sd_", bb,".prp.out")
outd = cbind(1:N, pvec)
write(file=outfile, t(outd), ncol=2)


