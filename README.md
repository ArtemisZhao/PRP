# Statistical Assessment of Replicability via Bayesian Model Criticism

## Description

This repository contains the software implementations for PRP. The software calculates the (prior/posterior)-replication predictive p-value, which is a statistical assessment of replicability via bayesian moodel criticism.

The repository includes source R code and scripts to replicate all the simulation results described in the manuscript. The data for real application is included in the R package.

## R source
Source code for R package `PRP` is included in `R_src`. To install, run

```{r}
devtools::install_github("ArtemisZhao/PRP/R_src")
```

## Simulation
The necessary code and data for simulation analysis are included in the folder __PRP_paper__. They should enable readers to fully reproduce the simulation results.

## Real data
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

## Contributors
Xiaoquan Wen (Umich)
Yi Zhao (Umich)

## References and Citations

