


# 
# 
# preprocess_fasta_database,
#  removes 1) duplicate entries (those with same sequence and accession)
# 		2) identical sequences if these come from the same species (retained if from different species)
#		3) seqeunces which are too long (note very long ones cause dereplication processes to crash)
#		4) sequences that are too short
# 
# 
# 
# 
# 
# 
# 
# 
# change log
#
# 24 Sep 2013: upper length limit option
# 10 Mar 2014: option to remove seqeunces with duplicate accessions.
# 11 Jun 2014: screen output for very large files.
# 09 Mar 2016: some settings now given by user in the command, instead of writing them in the script
# 12 May 2016: option -filter_all_duplicates, ignores test of whether from same species
# 14 Feb 2017: prints log of which are discarded
# 29 sep 2017: v10 of uclust
# 
# 
# 
# 
# 
# 
###########################################################################################################


my $arg_string  = join ' ', @ARGV;
#####################################
read_command_arguments($arg_string);#
#####################################




# Fung, Viri, Meta, Euka, Prok
# perl ~/usr_scripts/preprocess_fasta_database.pl
								# input:$database_file, output:$database_file.non_genome
# $remove_long_lines				= 0; 		# usearch dies if db contains sequences > about 25mb. good excuse to rm genome data
#	#$database_file				= "Hym.ng.rr.rep2.clusters.cl1";	# also remove v short seqs (< 40 bases)
#	$lowerlength_limit 			= 200;
#	$upper_length_limit 			= 32000;		# default:10000000
	$outfile_name 				= "$database_file.ng"; # non-genome


%printed_already;					
$duplicate_accessions_discarded;
								# uses usearch4.2.66_i86linux32
								# limited to about 2gb file.
								# input:$reduce_this_file, output:$reduce_this_file.rr
# $reduce_redundency				= 1;		# sequences discarded > 1mb. identical seqs removed, unless they come from different sp.
	$reduce_this_file			= "$database_file";
	$memory_efficient			= 1; 		# 1= read only one entry into memory at a time, allows large fasta databases to be used.


if($usearch_version == 4)
	{
	$dereplication_type 			= "derep_fullseq"; 
	};

if($usearch_version == 10)
	{
	$dereplication_type 			= "fastx_uniques"; 
	};

								# derep_subseq or derep_fullseq
								# sub seq uses a lot more memory..... but doesnt seem to work .. so ignore this option
								# it worked on hymeoptera (~100mb) but dies on insecta (~800mb)
	$remove_duplicate_accessions		= 1;					


$write_reduced_keyfile				= 0;		# first reads this file and stores the species identifier (string before underscore)
	$rewrite_keyfile_according_to_this_db	= "Euka.fas.ng.rr"; 
	$keyfile_to_rewrite			= "key_Aug2013_Euka";
	%species_level_tobycodes = ();
	$species_identifiers_used		= 3;		# 1 = tobycode, 2 = ncbi taxon number, 3= full species name



#################################################################################################################################
#
#
#
#
#
#################################################################################################################################




#################################################################################################################################
#
#
#
#
#
#################################################################################################################################





my %hits;
my $total_in_input = 0;
my $short_rm = 0;my $long_rm = 0;$not_rm=0;
my $bytes_rm =0;

if($remove_long_lines == 1)	{remove_long_lines();exit}
if($reduce_redundency == 1)	{reduce_redundency();exit}
if($write_reduced_keyfile == 1)	{write_reduced_keyfile();exit}










#################################################################################################################################
#
#
#
#
#
#################################################################################################################################



sub remove_long_lines
{
# usearch starts to cry if the database contains long sequences

print "\nreading db $database_file and removing seqs which are too short or too long\n";

open(FASTA_IN, $database_file) || die "Cant open $database_file.\n";
my $countentries_in_input=0;
while (my $line= <FASTA_IN>)
	{
	if($line =~ /^>/)
		{$countentries_in_input++}
	}
close(FASTA_IN);
print "
$database_file contains $countentries_in_input entries
";
open(FASTA_IN, $database_file) || die "Cant open $database_file.\n";
open(OUT1, ">$database_file.genome") || die "Cant open OUT.\n";
open(OUT2, ">$outfile_name") || die "Cant open OUT.\n";



my $fasta_entry = "";
my $countentries2=0;

while(my $fasta_line= <FASTA_IN>)
	{
#	print "$fasta_line";
	if($fasta_line =~ /^>.+/)
		{
		$countentries2++;
		if($countentries2 =~ /00000$/)
			{print "$countentries2 of $countentries_in_input, bytes removed:$bytes_rm\n"}

		unless(length($fasta_entry)<=2)
			{
			$fasta_entry =~ s/^>//;
			#print "\nfasta_entry:$fasta_entry\n";
			########################################
			process_entry($fasta_entry);#
			########################################
			}
		$fasta_entry = $fasta_line;
		}else{
		$fasta_entry .= $fasta_line;
		}
	};
close(FASTA_IN);

print "
total in input file:$total_in_input
removed because too short:$short_rm
removed because too long:$long_rm
printed to output:$not_rm
";


###############################################################################################

###############################################################################################



}





#################################################################################################################################
#
#
#
#
#
#################################################################################################################################




sub reduce_redundency
{
print "\n\nsub reduce_redundency\n";



print "\nsub reduce_redundency\n";


if($usearch_version == 4)
	{
# v4:
my $command = "$usearch_command --sort $reduce_this_file --output $reduce_this_file.sorted --maxlen 1000000";
system("rm $reduce_this_file.sorted");
print "command:$command\n";
system($command);

	}

# v10 seems not to use that



my $command;
if($usearch_version == 4)
	{
				#	--derep_subseq is more memory intensive than --derep_fullseq 
	# v 4
	$command = "$usearch_command --$dereplication_type --cluster $reduce_this_file.sorted --slots 4000003 --uc $reduce_this_file.uclust_log --seedsout $reduce_this_file.sorted.nr --maxlen 1000000";
	};

if($usearch_version == 10)
	{
	# v10:
	$command = "$usearch_command --$dereplication_type $reduce_this_file --slots 4000003 --uc $reduce_this_file.uclust_log";
	};


system("rm $reduce_this_file.uclust_log");
print "command:$command\n";
system($command);



my $identical_seq_same_sp =0;my $identical_seq_diff_sp =0;

open(IN, "$reduce_this_file.uclust_log") || die "\n\nerror 1779\n";
#open(OUT58, ">preprocess_fastadb_DEBUG");

open(OUTLOG , ">preprocess_LOG") || die "\nerror 235.\n";
print OUTLOG "Discarded\tIdential_To\n";

while (my $line = <IN>)
	{

	if($line =~ /^H\t.+\t100.0\t.+\t([\d\w].+)\t([\d\w].+)\n/)
		{
		my $hitseq = $1;my $hitseq2 = $2;my $tobcode = "unknown";my $tobcode2 = "unknown";
		if($hitseq =~ /(.+)_/){$tobcode = $1}
		if($hitseq2 =~ /(.+)_/){$tobcode2 = $1}
		if($tobcode eq $tobcode2)
			{
			$hits{$hitseq} = 1;$identical_seq_same_sp++;
			print OUTLOG "$hitseq\t$hitseq2\n";
			print OUT58 "identical sequence from same species to be rm:\n$line\n\n";
			}else{
			 $identical_seq_diff_sp++;#print "identical seq but different sp:$line\n\n";
			}
		if($filter_all_duplicates == 1)
			{
			$hits{$hitseq} = 1;
			};


		}

# H	358	15821	100.0	.	0	0	?	IDiB1T3T4D5D6Ba7ole_GU108467	IDiB1T3T4D5D6Ba7ole_GU108463

	}

close OUT58;
close OUTLOG;

my $count_discarded =0;
my $count_total =0;
my $count_printed_to_output =0;

my $file_as_string = "";
open(IN_FILTER, $reduce_this_file) || die "\nerror 1794 \n";
open(OUT, ">$reduce_this_file.rr") || die "\nerror 1327\n";

my $fasta_entry = "";
while (my $fasta_line = <IN_FILTER>)	
	{
	if($fasta_line =~ /^>.+/)
		{
		unless(length($fasta_entry)<=2)
			{
			$fasta_entry =~ s/^>//;


			########################################
			process_this_fasta_entry($fasta_entry);#
			########################################
			}
		$fasta_entry = $fasta_line;
		}else{
		$fasta_entry .= $fasta_line;
		}

	};

close(IN_FILTER);


#######################################################

sub process_this_fasta_entry
{
my $line = shift;
if($line =~ /^(.+)/)
	{
	my $speciesid = $1;	#print "$speciesid\n";
	if(exists($hits{$speciesid}) )
		{
		#print "discarding identical seq:$speciesid\n";
		$count_discarded++;
		}elsif(exists($printed_already{$speciesid}))
		{
		#print "discarding duplicate accession:$speciesid\n";
		$count_discarded++;
		$duplicate_accessions_discarded++;
		}else{
		$line =~ s/^.+\n//;
		print OUT ">$speciesid\n$line";$count_printed_to_output++;
		if($remove_duplicate_accessions ==1 ){$printed_already{$speciesid}=1};
		}
	$count_total++;
	}
}

#######################################################

print "\n\nentries in input ($count_total). discarded ($count_discarded), printed to output ($count_printed_to_output) \n";
print "duplicate accessions discarded:($duplicate_accessions_discarded)\n";
print "identical_seq_same_sp ($identical_seq_same_sp) identical_seq_diff_sp ($identical_seq_diff_sp)\n";


close(OUT);



$second_round = 0;

if($second_round == 1)
{
my $command = "$usearch_command --sort $reduce_this_file.rr --output $reduce_this_file.rr.sorted --maxlen 500000";
system("rm $reduce_this_file.rr.sorted");print "command:$command\n";system($command);
my $command = "$usearch_command --derep_subseq --cluster $reduce_this_file.rr.sorted --uc $reduce_this_file.rr.uclust_log --seedsout $reduce_this_file.rr.sorted.nr --maxlen 500000";
system("rm $reduce_this_file.rr.uclust_log");
print "command:$command\n";
system($command);


print "\nwarning, same-species-test might not be working on currently used format, needs checking.\n";


}







die;

open(IN8, "$reduce_this_file.rr") || die "\nerror 1828\n";

while(my $line = <IN8>)
	{
	if($line =~ />I([A-Z][a-z]{0,3})([A-Z][a-z]*.+)/)
#	if($line =~ />I(.+)/)
		{
#		my $tax = $1;
		my $order=$1; my $tax = $2;
	#	print "tax:$tax\n";
		$all_orders{$order}++;

		}
	if($line =~ />I([A-Z][a-z]{0,3})([A-Z][^0_4-9]+4)/)
		{
		my $familes=$2;
	#	print "familes:$familes\n";
		unless($line =~ /BOLD/){$all_familes{$familes}++};
		}

	if($line =~ />I([^0_7-9]+7)/)
		{
		my $genus=$1;
		unless($line =~ /BOLD/){$all_genus{$genus}++};
		}

	if($line =~ />I([^0_7-9]+7[^_]+)/)
		{
		my $sp=$1;
		unless($line =~ /BOLD/){$all_sp{$sp}++};
		}

	if($line =~ />I([^0_7-9]+7[a-z]+)/)
		{
		my $binom=$1;
		unless($line =~ /BOLD/){$all_binom{$binom}++};
		}

	if($line =~ />I([^_]*BOLD[^_]+)/)
		{
		my $bold=$1;
		$all_bold{$bold}++;
		}


	
	}

close(IN8);


my @orders  = keys %all_orders;
my @familes  = keys %all_familes;@familes = sort @familes;
my @genuss  = keys %all_genus;@genuss = sort @genuss;
my @sps  = keys %all_sp;@sps = sort @sps;
my @binoms  = keys %all_binom;@binoms = sort @binoms;
my @bolds  = keys %all_bold;@bolds = sort @bolds;


#foreach my $order(@orders){print "$order $all_orders{$order}\n"}
#foreach my $family(@familes){print "$family $all_familes{$family}\n"}
#foreach my $genus(@genuss){print "$genus $all_genus{$genus}\n"}
#foreach my $sp(@sps){print "$sp $all_sps{$sp}\n"}
#foreach my $binom(@binoms){print "$binom $all_binoms{$binom}\n"}
#foreach my $bold(@bolds){print "$bold $all_bold{$bold}\n"}


print "\n" , scalar @orders , " orders\n";
print "\n" , scalar @familes , " familes\n";
print "\n" , scalar @genuss , " genuss\n";
print "\n" , scalar @sps , " species\n";
print "\n" , scalar @binoms , " binom\n";
print "\n" , scalar @bolds , " bold\n";


###########################

# second round of derep, more memory intensive.... NOT USED


my $command = "usearch4.2.66_i86linux32 --sort $reduce_this_file.rr --output $reduce_this_file.rr.sorted --maxlen 500000";

#system("rm $reduce_this_file.rr.sorted");
#print "command:$command\n";
#system($command);
				#	--derep_subseq
#my $command = "usearch4.2.66_i86linux32 --derep_subseq --cluster $reduce_this_file.rr.sorted --uc $reduce_this_file.rr.uclust_log --seedsout $reduce_this_file.rr.sorted.nr --maxlen 500000";


#system("rm $reduce_this_file.rr.uclust_log");
#print "command:$command\n";
#system($command);








}






#################################################################################################################################
#
#
#
#
#
#################################################################################################################################




sub write_reduced_keyfile
{


###############
read_the_db();#
###############



#	$rewrite_keyfile_according_to_this_db	= "inv.fas.parsed.non_genome.rr";
#	$keyfile_to_rewrite			= "key_Mar2013_Insecta";


print "writing new key file ($keyfile_to_rewrite.reduced_according_to.$rewrite_keyfile_according_to_this_db) in which species that are not found in the database have been removed\n";


open (IN, "$keyfile_to_rewrite") || die "error 335. cant open file:$keyfile_to_rewrite\n";
open (OUT3, ">$keyfile_to_rewrite.reduced_according_to.$rewrite_keyfile_according_to_this_db") || die "error 336. cant open file\n";

my $count_keys = 0;
my $total_deleted =0;my $bold_deleted =0;
while ($line= <IN>)
	{
	# LZ1M4Mi7cal Micropterix_calthella 41027 species suborder:Zeugloptera family:Micropterigidae genus:Micropterix species:Micropterix calthella 

	if($line =~ /^(\S+)\s(\S+)\s(\d+)\sspecies /)
		{
		my $tobycode = $1;my $original_sp_id = $2;my $taxid = $3;

		if($species_identifiers_used == 1)# using tobycode as species label in fasta id's
			{
			if($species_level_tobycodes{$tobycode} == 2)
				{print OUT3 $line
				}else{
				$total_deleted++;if($original_sp_id =~ /BOLD/){$bold_deleted++}
				}
			}

			if($species_identifiers_used == 2)
				{# or using ncbi taxon number for species id's
				if($species_level_tobycodes{$taxid} == 2)
					{print OUT3 $line
					}else{
					$total_deleted++;if($original_sp_id =~ /BOLD/){$bold_deleted++}
					}
				}


		if($species_identifiers_used == 3)# using tobycode as species label in fasta id's
			{
			if($species_level_tobycodes{$original_sp_id} == 2)
				{print OUT3 $line
				}else{
				$total_deleted++;if($original_sp_id =~ /BOLD/){$bold_deleted++}
				}
			}



		}else{
		print OUT3 $line;
		}
	}


close(IN);
close(OUT3);

my $nonBOLDrm = $total_deleted-$bold_deleted;
print "total deleted ($total_deleted) bold deleted ($bold_deleted) nonBOLD removed ($nonBOLDrm) \n";

}



#################################################################################################################################
#
#
#
#
#
#################################################################################################################################



sub read_the_db
	{


#	$rewrite_keyfile_according_to_this_db	= "inv.fas.parsed.non_genome.rr";
#	$keyfile_to_rewrite			= "key_Mar2013_Insecta";



	print "sub fasta_format_to_nexus. reading $rewrite_keyfile_according_to_this_db. this may take a minute\n";

	my $count_missing_taxids =0;
	my $current_seq_length;
	my $missing_data_character = "N";

	open(FASTA_IN, $rewrite_keyfile_according_to_this_db) || die "Cant open input:$rewrite_keyfile_according_to_this_db\n";
	my $file_as_string = ""; my @all_lines = ();

	my $fasta_entry = "";
	while(my $fasta_line= <FASTA_IN>)
		{
		if ($fasta_line =~ /^>(.+)_([^_]+)\n/ )
			{
			$current_tobycode = $1;
			$species_level_tobycodes{$current_tobycode} = 2;
			}
		};
	close(FASTA_IN);


	};# 




#################################################################################################################################
#
#
#
#
#
#################################################################################################################################


sub process_entry
{
my $line = shift;
my $entrylength=length($line);
#print "\nprocessing:$line";
my $current_id = "";
if ($line =~ /^(.+)\n/ )
	{
	$current_id = $1;
	$line =~ s/^.+\n//;#$line =~ s/\n//;$line =~ s/\r//;my $current_seq_length = length($line);
	}else{
	die "\nerror 5631.\n"
	}

$total_in_input++;



if(length($line)<=$lowerlength_limit)
	{
#	print "why so short:($current_id) ($line)\n .... ignoring ... \n";
	$short_rm++;
	$bytes_rm+=$entrylength;

	}elsif(length($line)>=$upper_length_limit)
	{
#	print "why so long ($current_id):", length($line) , " .... ignoring ... \n";
	print OUT1 ">$current_id\n$line";
	$long_rm++;
	$bytes_rm+=$entrylength;
	}else{
	print OUT2 ">$current_id\n$line";
	$not_rm++;
	}

}



#################################################################################################################################
#
#
#
#
#
#################################################################################################################################




sub read_command_arguments
{

my $arguments = shift;
$arguments =~ s/\n//;$arguments =~ s/\r//;

my $error_message = "
\nsyntax error, quitting. \n
this script does 2 things, in order. Run like this:
perl preprocess_fasta_database.pl -in input_file_name -filter_seq_length_outliers -lower_length_limit 200 -upper_length_limit 32000
THEN
perl preprocess_fasta_database.pl -in input_file_name -reduce_redundency -usearch_command usearch4.2.66_i86linux32
\n";

print "
\n
";

if($arguments =~ /-in\s+(\S+)/)
	{
	$database_file = $1;
	}else{
	die "$error_message";
	};

if($arguments =~ /-filter_seq_length_outliers/)
	{
	$remove_long_lines	= 1;
	if($arguments =~ /-lower_length_limit\s+(\S+)/)
		{$lowerlength_limit = $1}else{die "$error_message"};
	if($arguments =~ /-upper_length_limit\s(\S+)/)
		{$upper_length_limit = $1}else{die "$error_message"};
	}elsif($arguments =~ /-reduce_redundency/)
	{
	$reduce_redundency = 1;	

	if($arguments =~ /-filter_all_duplicates/)
		{$filter_all_duplicates =1};

	if($arguments =~ /-usearch_command\s+(\S+)/)
		{$usearch_command = $1}else{die "$error_message"};


	my $test_usearch = `$usearch_command --help`;	# print "test_usearch:$test_usearch";
	unless($test_usearch =~ /Copyright 20.+Edgar/)
		{
		die "\nerror, no response from your usearch command ($usearch_command)\nis it installed and in path?\n"
		};

print "usearch_command:$usearch_command\n";

if($usearch_command =~ /usearch4\./)
	{
	$usearch_version = 4;print "\nappears you are using usearch version 4 ... good (has been tested a lot on this version)\n";
	}elsif($usearch_command =~ /usearch10\./)
	{
	$usearch_version = 10;print "\nappears you are using usearch version 10\n";
	}else{
	$usearch_version = 10;print "\nwarning, has only been tested on usearch v 4 and 10, could not detect either of these in your command...\n";
	
	};



# Copyright 2010-11 Robert C. Edgar



	#usearch4.2.66_i86linux32

	}else{
	die "$error_message";
	};

if($upper_length_limit >= 20000000){print "\nwarning, usearch does not accept very long seqeunces in DB, program might crash\n"};

# $remove_long_lines				= 0; 		# usearch dies if db contains sequences > about 25mb. good excuse to rm genome data
#	#$database_file				= "Hym.ng.rr.rep2.clusters.cl1";	# also remove v short seqs (< 40 bases)
#	$lowerlength_limit 			= 200;
#	$upper_length_limit 			= 32000;		# default:10000000

# $reduce_redundency

};


#################################################################################################################################
#
#
#
#
#
#################################################################################################################################


