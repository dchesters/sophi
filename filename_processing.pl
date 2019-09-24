
$in = $ARGV[0];


$TOL_version = 1;







open(IN, $in) || die "\nerror \n";

while (my $line   =<IN>)
	{
	$line =~ s/\n//;$line =~ s/\r//;
	if($line =~ /^(.+)/)
		{
		my $file  =$1;my $file2 = $file;my $file3 = $file;
	#	print "file2:$file2\n";
	#	die;

		if($TOL_version == 1)
			{
		# TOL version 1:
		$file2 =~ s/clo_pruned.fa.orth_ft_prune.+/clo_pruned.fa/;
		$command = "cp $file2 $file2.RM2";
			};

		if($TOL_version == 2)
			{
		# TOL version 2:
		# as named in file list:
		# phylobayesOUT.insectaNUCL.WP_100269.clo_sample.ppht
		# fasta file required:
		# insectaNUCL.WP_100066.clo_pruned.fa
		$file2 =~ s/phylobayesOUT.(insectaNUCL.WP\S+\.clo)_sample.ppht/$1_pruned.fa/;
		$command = "cp $file2 $file2.RMall";
			};


		if($TOL_version == 3)
			{
		# TOL version 3, using raxml gene trees:
	# in file list:
	# .....ections/nuclear_gene_trees/raxml_201606/RAxML_bipartitions.insectaNUCL.WP_413423.clo_pruned.fa.sed
	
	# fasta required:
	# insectaNUCL*.clo_pruned.fa


		$file2 =~ s/.+RAxML_bipartitions\.(insectaNUCL.+.clo_pruned.fa).sed/$1/;
		$command = "cp $file2 $file2.RM2";

			};




		if($TOL_version == 4)
			{
		$file2 =~ s/phylobayesOUT.(insectaNUCL.WP\S+\.clo)_sample.ppht/$1_pruned.fa.sed/;
		$command = "cp $file2 $file2.RMbin";
			};






		print "command:($command)\n";
		system($command);
		}

	}
close IN;






