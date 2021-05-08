all: simulation analysis additional plot	 
simulation: scripts/sim_batch_two_grp.R scripts/sim_batch_multi_grp.R
	if [ ! -d "sim_data" ]; then mkdir sim_data; fi
	perl scripts/sim_two_grp.pl
	perl scripts/sim_multi_grp.pl
analysis: scripts/analyze_batch_two_grp.R scripts/analyze_batch_multi_grp.R scripts/analyze_batch_multi_grp_metafor.R sim_data/ 
	if [ ! -d "output" ]; then mkdir output; fi
	perl scripts/analyze_two_grp.pl
	perl scripts/analyze_multi_grp_metafor.pl
	perl scripts/analyze_multi_grp.pl
additional: scripts/analyze_batch_multi_grp_fixeff.R scripts/sim_analyze_batch_two_grp_noisy.R
	Rscript scripts/sim_analyze_batch_two_grp_noisy.R 5
	Rscript	scripts/sim_analyze_batch_two_grp_noisy.R 10
	Rscript scripts/analyze_batch_multi_grp_fixeff.R 
plot: output/ 
clean:
	rm -rf sim_data/ output/
