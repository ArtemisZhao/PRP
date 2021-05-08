library(PRP)



args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
c_param = as.numeric(args[2])



d = read.table(input_file)

nc = dim(d)[2]

seq_b = seq(2,nc,2)
seq_s = seq(3,nc,2) 

beta = as.matrix(d[,seq_b])
se = as.matrix(d[,seq_s])


p = dim(d)[1]
pval_egger =c()
pval_Q = c()

for (i in 1:p){

	pval1 = posterior_prp(beta = beta[i,], se = se[i,], test=egger)$pval
	pval_egger = c(pval_egger,pval1)

	pval2 = posterior_prp(beta = beta[i,], se = se[i,], test=Q)$pval
	pval_Q = c(pval_Q, pval2)

}


filename = paste0("output/pubbias_multigrp_cparam_", round(c_param,1) , ".prp.out" )
outd = cbind(rep(1:p), pval_Q, pval_egger)
write(file=filename, t(outd), ncol=3)



