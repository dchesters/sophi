



# 
# 
# 
# 
# run job = 1 then job = 2, makes file of bash commands that you then run.
# 
# 
# 
# 
# 
# 
# change log
# 2017 JUL 11: seems ncbi changed directory structure for TSA's, required minor change here
# 	also option to not download TSA if its already in the working directory
# 2019 SEP 24: user specified variables transferred to command line
# 
# 
# 
# 
# 
# 
# 
# 
#############################################################################################





$in 		= $ARGV[0];
$job 		= $ARGV[1];
$limit_taxon 	= $ARGV[2];

unless($in =~ /[\w\d]/ && $job =~ /[12]/ && $limit_taxon =~ /\w+/)
	{
	die "\ncommand error.\n";
	};

# $job = 2; # 1 = download detail files. 2 = 






# $limit_taxon = "Hexapoda"; # Insecta




#####################################################################################################



if ($job == 1)
{
# read text format page of ftp://ftp.ncbi.nlm.nih.gov/genbank/tsa/
# find lines with tsa.+mstr.gbff.gz, these are links to summary on each project,
# download each summary 


open(IN, $in) || die "\nerror\n";
while (my $line = <IN>)
	{

	# ftp://ftp.ncbi.nlm.nih.gov/genbank/tsa/tsa.GAYY.mstr.gbff.gz

	if($line =~ /(tsa.+mstr.gbff.gz)/)
		{
		my $filename = $1;

# wont work any more:
#		my $wget_command = "wget -c ftp://ftp.ncbi.nlm.nih.gov/genbank/tsa/$filename";

		my $subdirectory;
		if($filename =~ /^tsa.(\w)/)
			{
			$subdirectory = $1;
			}else{
			die "\nerror cant parse first char to use as subdirectory\n"
			};

		my $wget_command = "wget -c ftp://ftp.ncbi.nlm.nih.gov/genbank/tsa/$subdirectory/$filename";


 		if (-e $filename) 
			{
	 		print "already have $filename\n";
 			}else{
			print "wget_command:$wget_command\n";
			system($wget_command);
			}

		};

	}
close IN;


};#if ($job == 1)




#####################################################################################################
#
#
#
#####################################################################################################




if ($job == 2)
{
print "\n\nJOB == 2, looking at summary files ... \n";

#/home/douglas/scripted_analyses/insect_TOL_analysis/data/tsa.GAAH.mstr.gbff.gz
 
my @file_list = split /\n/ , `ls tsa.G*.mstr.gbff*`;#.gz
@file_list = sort @file_list;
# print "\@file_list:@file_list\n";


open(OUT3, ">transcriptome_download_commands") || die "";


my $number_of_files = 0;my $number_of_insects = 0;my $store_command_index = 0;
foreach my $index (0 .. $#file_list)
	{
	my $file = $file_list[$index];
	my $unzipped_filename = $file;$unzipped_filename =~ s/\.gz$//;
#	print "\n$index of $#file_list, current file:$file unzipped_filename:$unzipped_filename\n";

#	print "\tcommand: (gunzip $file)\n";
	system("gunzip $file");

	my $current_file_is_insect = 0;my $found_genus = 0;my $genus = "";
	my $read_count_found = 0;my $read_count = 0;

	open(CURRENT, $unzipped_filename) || die "\nerror cant open ($unzipped_filename)\n";
	while (my $line2 = <CURRENT>)
		{
		$line2 =~ s/\n//;$line2 =~ s/\r//;
		# print "$line2\n";
		if($line2 =~ /$limit_taxon\;/)
			{$current_file_is_insect = 1};
		if($line2 =~ /  ORGANISM  ([A-Z][a-z]+) /)
			{$genus = $1;$found_genus =1};
		if($line2 =~ /LOCUS       \S+ +(\d+) rc +/)
			{$read_count = $1;$read_count_found = 1};

# number is count of reads
#LOCUS       GDQN01000000           11746 rc    RNA     linear   TSA 05-OCT-2016


		};
	close CURRENT;
#	print "\tprint running:($unzipped_filename)\n";
	system("gzip $unzipped_filename");


	if($current_file_is_insect == 1)
		{
		$number_of_insects++;
		print "\nFOUND $limit_taxon!\n";
		print "\t$index of $#file_list, current file:$file unzipped_filename:$unzipped_filename\n";


#current file:tsa.GAAB.mstr.gbff.gz unzipped_filename:tsa.GAAB.mstr.gbff
		my $seq_filename = $unzipped_filename;
		unless($seq_filename =~ s/^(tsa\.\w\w\w\w)\..+/$1/)
			{die "\neror 97, unzipped_filename:$unzipped_filename\n"};
		$seq_filename .= ".*.fsa_nt.gz";
		my $test_filename = $seq_filename; $test_filename =~ s/\*/1/;


		# Jul2017, directories have been changed:
		# ftp://ftp.ncbi.nlm.nih.gov/genbank/tsa/H/tsa.HAFV.1.fsa_nt.gz
		my $subdirectory;
		if($seq_filename =~ /^tsa.(\w)/)
			{$subdirectory = $1}else{die "\nerror cant parse first char to use as subdirectory\n"};


 		if (-e $test_filename) {
	 		print "already have $test_filename\n";
 			}else{
	 		print "dont have $test_filename\n";
			print OUT3 "echo $index of $#file_list\n";
		#	print OUT3 "wget ftp://ftp.ncbi.nlm.nih.gov/genbank/tsa/$seq_filename\n";
			print OUT3 "wget --tries=45 ftp://ftp.ncbi.nlm.nih.gov/genbank/tsa/$subdirectory/$seq_filename\n";
			};


	if($store_genus_names{$genus} == 1)
			{
			print "$genus AGAIN, read_count:$read_count, previous highest:$store_most_reads{$genus}\n";
			

			if($read_count =~ /\d+/ && $read_count >= $store_most_reads{$genus})
				{

				$genus_filtered_commands{$genus} = "echo $index of $#file_list\n" . 
					"wget --tries=45 ftp://ftp.ncbi.nlm.nih.gov/genbank/tsa/$subdirectory/$seq_filename\n";
				# hack:
			#	$genus_filtered_commands{$genus} = "cp $seq_filename ./genus_filtered/\n";


				};

			}else{
			$store_most_reads{$genus} = $read_count;

			$genus_filtered_commands{$genus} = "echo $index of $#file_list\n" . 
				"wget --tries=45 ftp://ftp.ncbi.nlm.nih.gov/genbank/tsa/$subdirectory/$seq_filename\n";
			# here hack since downloaded already, copy filtered to new folder:
		#	$genus_filtered_commands{$genus} = "cp $seq_filename ./genus_filtered/\n";


			};


	# ftp://ftp.ncbi.nlm.nih.gov/genbank/tsa/tsa.GAYY.mstr.gbff.gz

	#  2015 Aug 14        File        tsa.GCJV.mstr.gbff.gz  (1091 bytes)
	#  2015 Jun 11        File        tsa.GCJW.1.fsa_nt.gz  (19248006 bytes)

		$store_genus_names{$genus}=1;
		}
	$number_of_files++;

	#gzip file
	# file

	};

###################3


close OUT3;

	print "
	number_of_files:$number_of_files
	number_of_insects:$number_of_insects
	";

my @all_genera = keys %store_genus_names;
print "\ntotal number genera found:$#all_genera\n";


my @genus_filtered_commands_list = keys %genus_filtered_commands;
@genus_filtered_commands_list = sort @genus_filtered_commands_list;

open(OUT5, ">transcriptome_download_commands_genusfiltered") || die "";
foreach my $command(@genus_filtered_commands_list)
	{
	print OUT5 "$genus_filtered_commands{$command}";
	};
close OUT5;


};

#####################################################################################################
#
#
#
#####################################################################################################



#####################################################################################################
#
#
#
#####################################################################################################


