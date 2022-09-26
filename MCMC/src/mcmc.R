## random-walk MH: convergence issue difficult to cross barrier (0) if starting from the wrong side

# prior_PRP_mcmc <- function (hbo, so, hbr, sr, k=0.273, nreps = 50000, burnin_rate = 0.01, init=NA){
#
#
#    if(is.na(init)){
#        init = hbo
#    }
#    bbar = init
#    bbar_vec = c(bbar)
#
#
#
#    MH_step<-function( bbar, hbo, so, k){
#
#        # propose a new bbar via random walk proposal
#        bbar_p = rnorm(1, mean=bbar,sd=so)
#
#        #compute M-H ratio
#        # proposal ratio = 1
#        mhr = dnorm(hbo, mean = bbar_p, sd = sqrt(k^2*bbar_p^2+so^2))/dnorm(hbo, mean=bbar, sd=sqrt(k^2*bbar^2+so^2))
#        mhr = min(1, mhr)
#        if(rbinom(1,1,mhr)){
#            bbar = bbar_p
#        }
#
#        return(bbar)
#    }
#
#
#    for (i in 1:nreps){
#        bbar = MH_step(bbar,hbo,so,k)
#        bbar_vec = c(bbar_vec,bbar)
#    }
#
#    bbar_vec = bbar_vec[ceiling(nreps*burnin_rate):length(bbar_vec)]
#    p1 = sapply(bbar_vec, function(x) pnorm(hbr, mean=x, sd = sqrt(k^2*x^2+sr^2)))
#    p2 = 1 - p1
#    pval = 2*min(apply(cbind(p1,p2),2, mean))
#    cat("p-value = ", pval,"\n")
#    return(list(pval,bbar_vec))
#}



## independent MH algorithm

pmc_vec = c(0 , 0.001 , 0.002 , 0.003 , 0.004 , 0.005 , 0.006 , 0.007 , 0.008 , 0.009 , 0.01 , 0.011 , 0.012 , 0.013 , 0.014 , 0.015 , 0.016 , 0.017 , 0.018 , 0.019 , 0.02 , 0.021 , 0.022 , 0.023 , 0.024 , 0.025 , 0.026 , 0.027 , 0.028 , 0.029 , 0.03 , 0.031 , 0.032 , 0.033 , 0.034 , 0.035 , 0.036 , 0.037 , 0.038 , 0.039 , 0.04 , 0.041 , 0.042 , 0.043 , 0.044 , 0.045 , 0.046 , 0.047 , 0.048 , 0.049 , 0.05 )

k_vec = c(0 , 0.155368927298941 , 0.166005089694258 , 0.173337899205751 , 0.178978966221498 , 0.18388938695547 , 0.188067112764603 , 0.191791542513613 , 0.195271332191419 , 0.198575148690133 , 0.201534613519245 , 0.204481460061963 , 0.207059615021334 , 0.209648317091446 , 0.212212828958339 , 0.21443824547403 , 0.216731336214093 , 0.219075253455047 , 0.221090629773248 , 0.223007024495773 , 0.225129102734274 , 0.227105121873094 , 0.229073533229895 , 0.230927660776478 , 0.232684191110707 , 0.234490009720275 , 0.236500593128432 , 0.23807026226874 , 0.239798064432984 , 0.24137077435528 , 0.242982232444001 , 0.244687134655899 , 0.246358284045983 , 0.247862152521707 , 0.249551446935383 , 0.251016347829005 , 0.252623405238645 , 0.25419738927027 , 0.25548147869798 , 0.25715799531597 , 0.258514034625798 , 0.260067172562903 , 0.26149385602969 , 0.262818252318785 , 0.264347547831445 , 0.265740978439278 , 0.267262201879966 , 0.268573378227942 , 0.269859998040611 , 0.271219052445701 , 0.272612228503996)

k_max = max(k_vec)

prior_PRP_mcmc <- function (hbo, so, hbr, sr, nreps = 50000, burnin_rate = 0.01){


    bbar = hbo
    k = sample(k_vec,1)
    
    bbar_sample = c(bbar)
    k_sample = c(k)


    MH_step<-function( bbar, k, hbo, so ){

        # independent proposal
        k_p = sample(k_vec,1)
        bbar_p = rnorm(1, mean=hbo, sd = sqrt(k_max^2*hbo^2+so^2))


        #compute M-H ratio
        pr = dnorm(bbar_p, mean=hbo, sd = sqrt(k_max^2*hbo^2+so^2))/dnorm(bbar,  mean=hbo, sd = sqrt(k_max^2*hbo^2+so^2))
        lr = dnorm(hbo, mean = bbar_p, sd = sqrt(k_p^2*bbar_p^2+so^2))/dnorm(hbo, mean=bbar, sd=sqrt(k^2*bbar^2+so^2))
        mhr = lr/pr

        if(runif(1)<=mhr){
            bbar = bbar_p
            k = k_p
        }
        return(c(bbar,k))
    }


    for (i in 1:nreps){
        rst = MH_step(bbar,k,hbo,so)
        bbar_sample = c(bbar_sample,rst[1])
        k_sample = c(k_sample, rst[2])
        bbar = rst[1]
        k = rst[2]
    }

    bbar_sample = bbar_sample[ceiling(nreps*burnin_rate):length(bbar_sample)]
    p1 = sapply(1:length(bbar_sample), function(x) pnorm(hbr, mean=bbar_sample[x], sd = sqrt(bbar_sample[x]^2*k_sample[x]^2+sr^2)))
    p2 = 1 - p1
    pval = 2*min(apply(cbind(p1,p2),2, mean))
    cat("p-value = ", pval,"\n")
    
    return(pval)

}


