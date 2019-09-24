

@transcriptomes = @ARGV;

print "\@transcriptomes:@transcriptomes\n";





$path = "/home/douglas/scripted_analyses/insect_TOL_analysis/data/insect_TSA_processed2/";


foreach my $file(@transcriptomes)
	{
	$file_number++;
	print "file $file_number of $#transcriptomes, $file\n";

	my $peptide_file = "$path$file.transdecoder_dir/longest_orfs.pep";

	open(IN, $peptide_file) || die "\ncant open $peptide_file\n";

	my $file_as_string = "";
	while (my $line = <IN>)	{$file_as_string .= $line};
	close(IN);
	my @all_lines = split />/, $file_as_string;
	print "\t" , scalar @all_lines , " seqs in file\n";
	my %longest_orf = ();


	for my $each_line(1 .. $#all_lines)
		{
		my $line = $all_lines[$each_line];
		if($line =~ /^(.+)/)
			{
			my $speciesid = $1;	# print "$speciesid\n";
			$line =~ s/^.+\n//;
			$line =~ s/\012\015?|\015\012?//g;$line =~ s/\n//g;$line =~ s/\r//g;
			$line =~ s/[\s\t]//g;
			$line =~ s/\*$//;

			# Gene.5728::Agrilus_planipennis_429222167::g.5728::m.5728 type:complete len:113 gc:universal Agrilus_planipennis_429222167:166-504(+)

			my $sucinct_ID = "NA";
			if($speciesid =~ /^Gene\.\d+\:\:([^\:]+)/)
				{
				$sucinct_ID = $1
				}else{
				die "cant parse id from:$speciesid";
				};
		#	print "\tsucinct_ID:$sucinct_ID\n";

			my $orf_length = length($line);

			if(exists($longest_orf{$sucinct_ID}))
				{
				my $existing_length = length($longest_orf{$sucinct_ID});
			#	print "\tOLD ... existing_length:$existing_length orf_length:$orf_length\n";
				if($orf_length >= $existing_length)
					{
					$longest_orf{$sucinct_ID} = $line;
					}else{};
				}else{
				$longest_orf{$sucinct_ID} = $line;
			#	print "\tNEW ... \n";
				};

			};
		}



	my @keys = keys %longest_orf;

	open(OUT , ">$file.orf") || die "\nerro cant open out\n";

	foreach my $id(@keys)
		{	
		print OUT ">$id\n$longest_orf{$id}\n";
		unless($longest_orf{$id} =~ /^[A-Z]{10,}$/)
			{
			die "\nunexpected: file:$file id:$id longest_orf{id}:$longest_orf{$id}\n";
			};
		};

	close OUT;
	};









