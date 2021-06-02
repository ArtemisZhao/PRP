library(PRP)
library(metafor)
set.seed(123)

cat("COVID-19 Mortality Data\n\n")
data("mortality")

yi = mortality$beta
vi = mortality$se^2
rst=rma(yi,vi, method="FE")
erst = regtest(rst)

# posterior PRP 

# default
cat("Cochran's Q-test\n")
rst$QEp
cat("Posterior PRP: Q-quantity\n")
posterior_prp(beta=mortality$beta, se = mortality$se, test=Q, L=2000)$pvalue

cat("\n")
# publication bias
cat("Egger regression for detecting publication bias\n")
erst$pval
cat("Posterio PRP: Egger-quantitry\n")
posterior_prp(beta=mortality$beta, se = mortality$se, test=egger, L=2000)$pvalue
cat("\n")




cat("\nCOVID-19 Severity Data\n\n")
data("severity")
yi = severity$beta
vi = severity$se^2
rst=rma(yi,vi, method="FE")
erst = regtest(rst)



# default
cat("Cochran's Q-test\n")
rst$QEp
cat("Posterior PRP: Q-quantity\n")
posterior_prp(beta=severity$beta, se = severity$se, test=Q)$pvalue

cat("\n")

# publication bias
cat("Egger regression for detecting publication bias\n")
erst$pval
cat("Posterio PRP: Egger-quantitry\n")
posterior_prp(beta=severity$beta, se = severity$se, test=egger)$pvalue





