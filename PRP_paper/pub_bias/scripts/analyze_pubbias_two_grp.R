library(PRP)

args = commandArgs(trailingOnly=TRUE)
input_file = args[1]
p_thresh = as.numeric(args[2])



d = read.table(input_file)
beta = as.matrix(cbind(d$V2, d$V4))
se = as.matrix(cbind(d$V3,d$V5))
p = length(d$V2)

pval1 = sapply(1:p, function(x) prior_prp(beta = beta[x,], se = se[x,]),test="pub_bias")
pval2 = sapply(1:p, function(x) prior_prp(beta = beta[x,], se = se[x,]),test="two_sided")
filename = paste0("output/pubbias_2grp_pthresh_", round(p_thresh,2) , ".prp.out" )
outd = cbind(rep(1:p), pval2, pval1)
write(file=filename, t(outd), ncol=3)

