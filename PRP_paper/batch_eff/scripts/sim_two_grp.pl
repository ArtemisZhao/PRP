@bb = (0, 0.2, 0.4, 0.6, 0.8, 1.0);

foreach $b (@bb){

	print "Rscript scripts/sim_batch_two_grp.R  $b\n";
	`Rscript scripts/sim_batch_two_grp.R  $b`;
}

