library(PRP)


args = commandArgs(trailingOnly=TRUE)
filename = args[1]

bb = args[2]




d = read.table(filename)
attach(d)
N = dim(d)[1]

pvec = sapply(1:N, function(x) prior_prp(beta=c(V2[x], V4[x]), se=c(V3[x], V5[x]), test="two_sided")$pvalue)

outfile = paste0("output/batch_2grp_bb_sd_",bb,".prp.out")
outd = cbind(1:N, pvec)
write(file=outfile, t(outd), ncol=2)


