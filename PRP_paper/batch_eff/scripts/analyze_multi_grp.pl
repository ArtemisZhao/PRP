@files = <sim_data/batch_multigrp*.dat>;
foreach $f (@files){

	next if $f !~ /sd\_(\S+)\.dat/;
	$b = $1;
	print "Rscript scripts/analyze_batch_multi_grp.R  $f $b\n";
	`Rscript scripts/analyze_batch_multi_grp.R  $f $b`;
}

