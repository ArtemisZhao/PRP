all: simulation analysis plot	 
simulation: scripts/sim_pubbias_two_grp.R scripts/sim_pubbias_multi_grp.R
	if [ ! -d "sim_data" ]; then mkdir sim_data; fi
	Rscript scripts/sim_pubbias_two_grp.R 1
	Rscript scripts/sim_pubbias_two_grp.R 0.05
	Rscript scripts/sim_pubbias_two_grp.R 0.01
	Rscript scripts/sim_pubbias_multi_grp.R  0
	Rscript scripts/sim_pubbias_multi_grp.R  5
	Rscript scripts/sim_pubbias_multi_grp.R  10
analysis: scripts/analyze_pubbias_two_grp.R scripts/analyze_pubbias_multi_grp.R scripts/analyze_pubbias_multi_grp_metafor.R sim_data/ 
	if [ ! -d "output" ]; then mkdir output; fi
	Rscript scripts/analyze_pubbias_two_grp.R sim_data/pubbias_2grp_pthresh_1.dat 1
	Rscript scripts/analyze_pubbias_two_grp.R sim_data/pubbias_2grp_pthresh_0.05.dat 0.05
	Rscript scripts/analyze_pubbias_two_grp.R sim_data/pubbias_2grp_pthresh_0.01.dat 0.01
	Rscript scripts/analyze_pubbias_multi_grp_metafor.R sim_data/pubbias_multigrp_cparam_0.dat 0
	Rscript scripts/analyze_pubbias_multi_grp_metafor.R sim_data/pubbias_multigrp_cparam_5.dat 5
	Rscript scripts/analyze_pubbias_multi_grp_metafor.R sim_data/pubbias_multigrp_cparam_10.dat 10
	Rscript scripts/analyze_pubbias_multi_grp.R sim_data/pubbias_multigrp_cparam_0.dat 0 &
	Rscript scripts/analyze_pubbias_multi_grp.R sim_data/pubbias_multigrp_cparam_5.dat 5 &
	Rscript scripts/analyze_pubbias_multi_grp.R sim_data/pubbias_multigrp_cparam_10.dat 10 
plot: output/ 
clean:
	rm -rf sim_data/ output/

