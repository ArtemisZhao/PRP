## Simulation analysis

This directory contains the code to reproduce the simulation results presented in the paper. 
Each sub-directory represents a set of analyses (the batch effect simulation and the publication bias simulation). 
To generate the results, simply go into each directory, and run

```{r}
make
```

This will start the process of 
data generation -> data analysis -> summarizing results. 


In all cases, the intermediate results (simulated data) are saved in the generated data directory. 
The final summary of results are saved in the generated results directory. 

To restore the original directory structure before analysis, run

```{r}
make clean
```
