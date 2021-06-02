# Statistical Assessment of Replicability via Bayesian Model Criticism

## Description

This repository contains the software implementations for PRP. The software calculates the (prior/posterior)-replication predictive p-value, which is a statistical assessment of replicability via bayesian moodel criticism.

The repository includes source R code and scripts to replicate all the simulation results described in the manuscript. The data for real application is included in the R package.

## R source
Source code for R package `PRP` is included in `R_src`. To install, run

```{r}
devtools::install_github("ArtemisZhao/PRP/R_src")
```

## Data and code
The code and data for simulation and real data analysis are included in the folder __PRP_paper__. They should enable readers to fully reproduce the simulation results.

To reproduce simulation and analysis results, clone the repo and run ```make``` in each sub-directory. The simulation data will be re-generated. 

The data for the real application described in the paper are included in the R package.

```{r}
library(PRP)
## RP:P 
data("RPP_filtered")

## Cardiovascular disease impact on the COVID-19 mortality
data("mortality")

## Cardiovascular disease impact on the COVID-19 severity
data("severity")
```

### R library dependency

The following R packages are required to install and run the analysis code:

+ ```devtools```
+ ```mvtnorm```
+ ```metafor```
+ ```dplyr```
+ ```ggplot2```

## Docker image

A docker image with pre-configured Linux running environment and pre-installed R libraries is availabel for download from [docker hub](https://hub.docker.com/r/xqwen/prp).


## Contributors

+ Yi Zhao (zhayi at umich dot edu)
+ Xiaoquan Wen (xwen at umich dot edu)

## References and Citations

+ Zhao, Y. and Wen, X. Statistical Assessment of Replicability via Bayesian Model Criticism. [arXiv:2105.03993](https://arxiv.org/abs/2105.03993) 
