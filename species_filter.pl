
#
#
#
#
#	species_filter.pl 
#	perl script for filtering fasta file according to species names (optimized for idosyncrasies of insect labels)
#
#
#    	Copyright (C) 2017  Douglas Chesters
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#	contact address dc0357548934@live.co.uk
#
#
#
#
#
# 
# #############################################################################################
# 
#
#
#	options:
# 			-in -format -trim_accession -identified_species_only -filter_method
#
#
#
#
#
# 
# 
# change log 
# 01Aug2014: for each species, select the most representative (as opposed to the longest)
# 04Aug2014: bugfix for printing identified / non-identified
# 05Aug2014: option to filter by least ambiguous bases.
# 26Dec2014: insignificant changes
# 16Jan2015: deletes subspecies names
# 24Jan2015: discard entries where cant parse binomial.
#		removes substrings for verbose species IDs
# 31Jan2015: stop printing 'nr' species 
# 10Apr2016: options given in command, not written in script.
# 10Aug2017: since ncbi discontinued GI numbers, made modifications for using accessions
# 		other corrections to regexes, to stop OTU labels being overfiltered
# 03Nov2017 (v1.02): 
#		completly inactivated regex's for parsing GI numbers
# 		newline read bugfix 
# 		reports counts for different regexes used in parsing
# 		added GPL (curious it didnt have one already, 
#			surely this would have been in at least one of the published pipelines)
# 10Sep2019: an option for strict binomial filtering.
# 
#
# 
# 
# 
##################################################################################################



my $arg_string  = join ' ', @ARGV;
#####################################
read_command_arguments($arg_string);#
#####################################


my $version = 1.02;


# my $input 			= $ARGV[0];
# my $format 			= $ARGV[1]; #if($format == 2)species_filtering_tobycoded else species_filtering();

# $trim_off_accession 		= 0;	# must be 0 or 1

# $identified_species_only 	= 1;
# $filter_by_MRS			= 0; 	# 0 = take longest sequence for each species. 
					# 1 = take the most representative sequence for each species.
					# 2 = take the sequence with least ambiguous data (N's etc). these make alignment more difficult. 

$MRS_script 			= "~/usr_scripts/most_representative_seq.pl";



$print_screen			= 1;	# 0 = quiet

############################################################################################################



$date = localtime time;

open(LOG, ">>species_filter_LOG");
print LOG "\n$date input:$input format:$format\n";
system("rm $input.ID_filtered");


unless($print_screen == 0)
{
print "
\n ***** species_filter.pl (v $version)  *****  
";
if($filter_by_MRS == 0){print "filter_by_MRS == 0. take longest sequence for each species.\n"};
if($filter_by_MRS == 1)
	{
	print "filter_by_MRS == 1. take the most representative sequence for each species.\n";
	print "MRS_script:$MRS_script\n";
	};
if($filter_by_MRS == 2){print "filter_by_MRS == 2. take the sequence with least ambiguous data (N's etc). these make alignment more difficult.\n"};
print "trim_off_accession:$trim_off_accession\n";
print "identified_species_only:$identified_species_only\n\n";

};



if($no_accession == 1)
	{
	filter_no_accessions();
	die "";
	};





if($format == 2)
	{
	species_filtering_tobycoded();
	}else{
	species_filtering();
	}

close(LOG);

if($regex_hits{1} =~ /\d/){print "regex 1 (^[A-Z][a-z]+_sp_[^\|\_]+_[A-Z]{2}[0-9]{6}),\n\tinstances found:$regex_hits{1}\n"};
if($regex_hits{2} =~ /\d/){print "regex 2 (^[A-Z][a-z]+_sp_BOLD_[\w\d]+_[A-Z]{2}[0-9]{6}+),\n\t instances found:$regex_hits{2}\n"};
if($regex_hits{3} =~ /\d/){print "regex 3 (^[A-Z][a-z]+_sp_[a-zA-Z0-9\_]+_[A-Z]{2}[0-9]{6}+),\n\t instances found:$regex_hits{3}\n"};
if($regex_hits{4} =~ /\d/){print "regex 4 (^[A-Z][a-z]+_aff_[^\|\_]+_[A-Z]{2}[0-9]{6}+),\n\t instances found:$regex_hits{4}\n"};
if($regex_hits{5} =~ /\d/){print "regex 5 (^[A-Z][a-z]+_[a-z]+_[a-z]+_[A-Z]{2}[0-9]{6}+),\n\t instances found:$regex_hits{5}\n"};
if($regex_hits{6} =~ /\d/){print "regex 6 (^[A-Z][a-z]+_[a-z]+_[A-Z]{2}[0-9]{6}+),\n\t instances found:$regex_hits{6}\n"};
if($regex_hits{7} =~ /\d/){print "regex 7 (^[A-Z][a-z]+_[A-Z][a-z]+_sp_.+_[A-Z]{2}[0-9]{6}),\n\t instances found:$regex_hits{7}\n"};
if($regex_hits{8} =~ /\d/){print "regex 8 (^.+[\|\_][^\|\_]+),\n\t instances found:$regex_hits{8}\n"};












unless($print_screen == 0){print "\nFIN.\n"};


exit;

##############################################################################################################
#
#
#
#
#
##############################################################################################################




sub species_filtering
{

unless($print_screen == 0){print "\nsub species_filtering\n\n"};

$discarded_entries = 0;

open (IN_FILTER, $input) || die "cant open $input\n";
my $file_as_string = "";
while (my $line = <IN_FILTER>){$file_as_string .= $line}
close(IN_FILTER);

my @all_lines = split />/, $file_as_string;
unless($print_screen == 0){print scalar @all_lines , "\tseqs in your input file.\n"};



for my $each_line(1 .. $#all_lines)
	{
	my $line = $all_lines[$each_line];#	if($line =~ /^([^\_]+[_\.][^\_]+)(.*)/)
	if($line =~ /(.+)/)
		{
		my $description = $1;$description =~ s/[\.\-]//g;#print "\ndescription:$description\n";
		$description =~ s/\n//;$description =~ s/\r//;
		my $speciesid = "NA";my $rest = "NA";
		my $discard_entry = 0;
		my $regex_hit = "NA";


		# following for accessions
		if($description =~ /^([A-Z][a-z]+_sp_[^\|\_]+)(_[A-Z]{2}[0-9]{6})$/)
			{

			$speciesid = $1;$rest = $2;$regex_hit=1;$regex_hits{$regex_hit}++;# print "speciesid:$speciesid accession:$rest\n";
			if($identified_species_only == 1){$discard_entry = 1;$discarded_entries++};

			# >Collembola_sp_BOLD_ACD2633_KM619749
			}elsif($description =~ /^([A-Z][a-z]+_sp_BOLD_[\w\d]+)(_[A-Z]{2}[0-9]{6}+)$/){
 				 $speciesid = $1;$rest = $2;$regex_hit = 2;$regex_hits{$regex_hit}++;
				if($identified_species_only == 1){$discard_entry = 1;$discarded_entries++};

			# >Ablabesmyia_sp_20225_C1_KF000166
			}elsif($description =~ /^([A-Z][a-z]+_sp_[a-zA-Z0-9\_]+)(_[A-Z]{2}[0-9]{6}+)$/){
 				 $speciesid = $1;$rest = $2;$regex_hit = 3;$regex_hits{$regex_hit}++;
				if($identified_species_only == 1){$discard_entry = 1;$discarded_entries++};

			}elsif($description =~ /^([A-Z][a-z]+_aff_[^\|\_]+)(_[A-Z]{2}[0-9]{6}+)$/){
 				 $speciesid = $1;$rest = $2;$regex_hit = 4;$regex_hits{$regex_hit}++;
				if($identified_species_only == 1){$discard_entry = 1;$discarded_entries++};

			}elsif($description =~ /^([A-Z][a-z]+_[a-z]+)_([a-z]+)(_[A-Z]{2}[0-9]{6}+)$/)#Apis_mellifera_intermissa_693583962
			{
			$speciesid = $1; my $subspecies = $2; my $gi=$3;$rest = $gi;$regex_hit = 5;$regex_hits{$regex_hit}++;

			}elsif($description =~ /^([A-Z][a-z]+_[a-z]+)(_[A-Z]{2}[0-9]{6}+)$/) # species_ID:(Callytron_yuasai) rest:(_KC963819)
			{
			$speciesid = $1;$rest = $2;$regex_hit = 6;$regex_hits{$regex_hit}++;
			}elsif($description =~ /^([A-Z][a-z]+_[A-Z][a-z]+_sp_.+)(_[A-Z]{2}[0-9]{6}$)/)
			{
			$speciesid = $1;$rest = $2;$regex_hit = 7;$regex_hits{$regex_hit}++;
			# >Cicindela_Rivacindela_sp_JP_102_AJ618260


	#	# following for GI numbers:
	#	}elsif($description =~ /^([A-Z][a-z]+_sp_[^\|\_]+)(_\d+)$/)#
	#		{
	#		$speciesid = $1;$rest = $2;$regex_hit = 7;$regex_hits{$regex_hit}++;
	#		if($identified_species_only == 1){$discard_entry = 1;$discarded_entries++};
	#		}elsif($description =~ /^([A-Z][a-z]+_aff_[^\|\_]+)(_\d+)$/){
 	#			 $speciesid = $1;$rest = $2;$regex_hit = 8;$regex_hits{$regex_hit}++;
	#			if($identified_species_only == 1){$discard_entry = 1;$discarded_entries++};
	#		}elsif($description =~ /^([A-Z][a-z]+_[a-z]+)_([a-z]+)(_\d+)$/)#Apis_mellifera_intermissa_693583962
	#		{
	#		$speciesid = $1; my $subspecies = $2; my $gi=$3;$rest = $gi;$regex_hit = 9;$regex_hits{$regex_hit}++;
	#		#print "speciesid:$speciesid subspecies:$subspecies gi:$gi\n";
	#		}elsif($description =~ /^([A-Z][a-z]+_[a-z]+)(_.+)$/)#Helicoverpa_punctigera_586947476
	#			{
	#			$speciesid = $1;$rest = $2;$regex_hit = 10;
	#			if($rest =~ s/.+_.+_.+_(.+)$/$1/)# must be very long id
	#				{
	#			#	print "\nwarning. speciesid:$speciesid, long ID ($description) doesnt all look like accession, ";
	#			#	print "will cut some of the string out, retaining just this:($rest) \n"
	#				};
	#			#>Colletes_fasciatus_s_l_NML_2007_EF028541
	#			#>Colletes_clypearis_group_sp_ZQN_2013_KC469666

				}else{
				if($identified_species_only == 1)
					{
					print "\nWARNING 174. you have specified binomials only, and I cant parse binomial from ID:($description) ... SO discarding this entry.\n";
					$discard_entry = 1;$discarded_entries++;
					}else{
				#	print "\nWARNING 177. I cant parse binomial from ID:($description) ... ";
					# Cicindela_Calochroa_cf_safraneki_KT_2015_LC020304
				#	print "although you have specified retain non-binomials (\$identified_species_only = $identified_species_only).\n";
					if($description =~ /^(.+)([\|\_][^\|\_]+)$/)
						{
						$speciesid = $1;$rest = $2;	$regex_hit = 8;$regex_hits{$regex_hit}++;

						}else{die "\ni dont know what this is ($description). quitting\n"};
					};

				};



		
	#	print "\nregex_hit:($regex_hit) species_ID:($speciesid) rest:($rest)\n";
		$line =~ s/^.+\n//;
		$line =~ s/\n//g;$line =~ s/\r//g;$line =~ s/\012\015?|\015\012?//g;

		unless($discard_entry == 1)
		{

		unless($all_accessions_this_species{$speciesid} =~ /$rest\t/)
			{
			$all_accessions_this_species{$speciesid} .= "$rest\t"; 
			};


		if($filter_by_MRS == 1)
			{

			$sequences{$speciesid} .= ">$speciesid$rest\n$line\n";

			}else{



			if($filter_by_MRS == 0)	# take longest
				{
				if(exists($sequences{$speciesid}))
					{
					my $len_existing = length($sequences{$speciesid});
					my $len_new = length($line);	#	print "\texisting:$len_existing new:$len_new\n";
					if($len_new >= $len_existing)
						{$sequences{$speciesid} = $line;$therest{$speciesid} = $rest;#print "\treplacing\n"
						}
					}else{
					$sequences{$speciesid} = $line;
					$therest{$speciesid} = $rest;
					}
				}


			if($filter_by_MRS == 2)	# take that with least ambig data
				{
				if(exists($sequences{$speciesid}))
					{
					print LOG "speciesid:$speciesid observed already\n";

					my $existing_seq = $sequences{$speciesid};my $new_seq = $line;
					#print "existing_seq:$existing_seq\nnew_seq:$new_seq\n";
					my $count_ambig_existing = 0;my $count_ambig_new = 0;
					while($existing_seq =~ s/[nyr]//i){$count_ambig_existing++};
					while($new_seq =~ s/[nyr]//i){$count_ambig_new++};
					print LOG "\tcount_ambig_existing:$count_ambig_existing count_ambig_new:$count_ambig_new\n";

					if($count_ambig_new < $count_ambig_existing)
						{$sequences{$speciesid} = $line;$therest{$speciesid} = $rest;#print "\tless ambiguous data. replacing\n"
						}

					}else{
					print LOG "speciesid:$speciesid is novel\n";
					$sequences{$speciesid} = $line;
					$therest{$speciesid} = $rest;
					
					};				
				}


			}

			};

		}else{
		print "\n\nwarning. unexpected ID.\n"
		}
	}

my @unique_species = keys %sequences;
@unique_species = sort @unique_species;
open(OUT, ">$input.ID_filtered") || die "\n\nerror\n";
open(OUT2, ">$input.retained_accessions") || die "\n\nerror\n";
close(OUT);


my $printed_to_output =0;


foreach my $sp(@unique_species)
	{


	if($filter_by_MRS == 1)
		{
		

		my $theMRS = find_MRS($sp);

			open(OUT, ">>$input.ID_filtered") || die "\n\nerror\n";
			if($trim_off_accession == 1){$theMRS =~ s/(>[A-Z][a-z]+_[a-z]+).+/$1/}
			print OUT "$theMRS";
			close OUT;
$printed_to_output++;

		}else{

		my $print_current = "";
		if($identified_species_only == 1)
			{
			if($sp =~ /^[A-Z][a-z]+[_][a-z][a-z]+$/ )#>Ceratina_n_subgen_sp_SMR_2010_GU321539
				{
				# >Euglossa_cf_variabilis_SR417_AY920314

				if( $sp =~ /^[A-Z][a-z]+[_]sp$/ || $sp =~ /^[A-Z][a-z]+[_]aff$/ || 
					$sp =~ /^[A-Z][a-z]+[_]cf$/ || $sp =~ /^[A-Z][a-z]+[_]nr$/)
					{
					$print_current = 0;
					}else{
					$print_current = 1;
		
					}
				
				#if($sp =~ /Braunsapis_nr/){print "sp:$sp\n";die "\n\n"}

				}else{
				my $print_current = 0;
				}  
			}else{
			$print_current =1;
			}

		if($print_current == 1)
			{
		#	print "sp:$sp\n";
			
			open(OUT, ">>$input.ID_filtered") || die "\n\nerror\n";
			print OUT ">$sp" , "";
			if(length($therest{$sp})>=2 && $trim_off_accession == 0)
				{print OUT $therest{$sp}}
			print OUT "\n" , $sequences{$sp}  , "\n";#	print ">$sp\n" , $sequences{$sp} , $therest{$sp} , "\n";
			close OUT;
			$printed_to_output++;

			my $all_accessions = $all_accessions_this_species{$sp};
			print OUT2 "$sp\t$therest{$sp}\t$all_accessions\n";

			};

		}

	}#foreach my $sp(@unique_species)

close OUT2;

unless($print_screen == 0)
	{
print "" , scalar @unique_species , "\tunique_species,
$discarded_entries\tdiscarded_entries (only relevent if user specified -identified_species_only),
$printed_to_output\tprinted_to_output.\n
"
	};


}



############################################################
#
#
#
#
#
############################################################


sub species_filtering_tobycoded
	{

print "\nspecies_filtering, assuming taxstring underscore accession\n";


open (IN_FILTER, $input) || die "cant open $input\n";

my $file_as_string = "";

while (my $line = <IN_FILTER>)
	{$file_as_string .= $line}
close(IN_FILTER);

my @all_lines = split />/, $file_as_string;

print scalar @all_lines , " seqs in file\n";

for my $each_line(1 .. $#all_lines)
	{
	my $line = $all_lines[$each_line];

	if($line =~ /^(.+)[_]([^\_\n\r]+)/)
		{
		my $speciesid = $1;my $rest = $2;$rest =~ s/[\n\r \s]//g;unless($rest =~ /^[A-Z]+[0-9]+$/){die "\n\nwierd accession in:$rest\n"}
	#	print "$speciesid rest:$rest\n";
		$line =~ s/^.+\n//;
		$line =~ s/\n//g;$line =~ s/\r//g;$line =~ s/\012\015?|\015\012?//g;

		if(exists($sequences{$speciesid}))
			{
			my $len_existing = length($sequences{$speciesid});
			my $len_new = length($line);
		#	print "\texisting:$len_existing new:$len_new\n";
			$number_replacments{$speciesid}++;

			if($len_new >= $len_existing){$sequences{$speciesid} = $line;$therest{$speciesid} = $rest;
				#print "\treplacing\n"
				}
			}else{
			$sequences{$speciesid} = $line;
			$therest{$speciesid} = $rest;
			}

		}
	}

my @unique_species = keys %sequences;
@unique_species = sort @unique_species;
open(OUT, ">$input.ID_filtered") || die "\n\nerror\n";

my $countprinted=0;

foreach my $sp(@unique_species)
	{

	my $print_current = 1;
	if($identified_species_only == 1)
		{
		unless($sp =~ /^[A-Z][a-z]+[_][a-z]+$/){$print_current = 0}# only works with format 1
		}else{
		$print_current =1;
		}


	if($print_current == 1)
		{
		print OUT ">$sp" , "";
	$countprinted++;
		if(length($therest{$sp})>=2)
			{
			print OUT $therest{$sp}; 
			}

		print OUT "\n" , $sequences{$sp}  , "\n";
	#	print ">$sp\n" , $sequences{$sp} , $therest{$sp} , "\n";
		}

if($number_replacments{$sp}>=100){print LOG "species ID:$sp occured:$number_replacments{$sp} times in file\n"}

unless($sp =~ /^[^7]+7[a-z]+$/){$unidentified++}
	}


close(OUT);

print scalar @all_lines , " seqs in file, filtered to:$countprinted\n";
print LOG scalar @all_lines , " seqs in file, filtered to:$countprinted. unidentified:$unidentified\n";



}



############################################################
#
#
#
#
#
############################################################


sub find_MRS
{
my $species = shift;

my $data = $sequences{$species};


my @all_lines = split />/, $data;splice @all_lines , 0,1;
unless($print_screen == 0){print "species:$species, seqs:" , scalar @all_lines , "\n"};

my $MRS = "NA";

if(scalar @all_lines == 1)
	{
	
	$MRS = $data;
#	print "only 1 sequence:$MRS\n";

	}elsif(scalar @all_lines == 2)
	{
	# 2 sequences for current species, pick at random.

	$MRS = ">$all_lines[0]";
#	print "only 2 sequence, selecting 1 randomly:$MRS\n";

#	die;
	}else{

	# 3 or more sequences, find MRS

	open(OUT77, ">all_seqs_current_sp") || die "\nerror 286\n";
	print OUT77 "$data";
	close OUT77;

	my $command = "perl $MRS_script -i all_seqs_current_sp -s 120 -b blastn-2.2.28+-64bit -a 1200 -p 80";
	# -i followed by name of fasta file of DNA sequences
	# -s followed by number of sequences to sample from the input file,
	#    from which the MRS will be found
	# -b followed by the name of the command for running blast on your system
	#    e.g. blastn or /home/plantID/ncbi-blast-2.2.27+/bin/blastn or blastn-2.2.28+-64bit
	# -a followed by number, which is sequence length below which hits are ignored
	# -p follows by number, percent identity cutoff

	system("rm all_seqs_current_sp.MRS");
	system($command);

	my $test_MRS = `cat all_seqs_current_sp.MRS`;

	print "test_MRS:($test_MRS)\n";
	my @test_splitfile = split />/, $test_MRS;splice @test_splitfile , 0,1;

	if(scalar @test_splitfile == 1)
		{
		$MRS = 	$test_MRS;		
		}else{
		print "\nerror:MRS output contains more than one sequence. or none. just taking first by random.\n";
		$MRS = ">$all_lines[0]";

		}
	#die;
	}

return($MRS);

}






############################################################
#
#
#
#
#
############################################################



##############################################################################################################




sub read_command_arguments
{

my $arguments = shift;
$arguments =~ s/\n//;$arguments =~ s/\r//;

print "
";


if($arguments =~ /-in\s+(\S+)/)
	{
	$input = $1;
	}else{
	die "\n\nerror, please give input \n";
	};
# my $input 			= $ARGV[0];

if($arguments =~ /-format\s+(\d)/)
	{
	$format = $1;
	}else{
	die "\n\nerror 673\n";
	};
# my $format 			= $ARGV[1]; #if($format == 2)species_filtering_tobycoded else species_filtering();



if($arguments =~ /-trim_accession/)
	{
	$trim_off_accession = 1;
	}else{
	$trim_off_accession = 0;
	};
# $trim_off_accession 		= 0;	# must be 0 or 1


if($arguments =~ /-no_accession/)
	{
	$no_accession = 1;
	};



if($arguments =~ /-identified_species_only/)
	{
	$identified_species_only = 1;
	}else{
	$identified_species_only = 0;
	};
# $identified_species_only 	= 1;

if($arguments =~ /-filter_method\s+(\d)/)
	{
	$filter_by_MRS = $1;
	}else{
	die "\n\nerror 569\n";
	};
# $filter_by_MRS			= 0; 	# 0 = take longest sequence for each species. 
					# 1 = take the most representative sequence for each species.
					# 2 = take the sequence with least ambiguous data (N's etc). these make alignment more difficult. 



# -in -format -identified_species_only -filter_method 2




};

##############################################################################################################




sub filter_no_accessions
{


	open(FASTA_IN, $input) || die "Cant open $database_file.\n";
	open(OUT1, ">$input.binom_only") || die"";
	my $fasta_entry = "";
	while(my $fasta_line= <FASTA_IN>)
		{
		if($fasta_line =~ /^>.+/)
			{
			unless(length($fasta_entry)<=2)
				{
				$fasta_entry =~ s/^>//;
				########################################
				process_entryF($fasta_entry);#
				########################################
				}
			$fasta_entry = $fasta_line;
			}else{
			$fasta_entry .= $fasta_line;
			}
		};
	close(FASTA_IN);

	# do the last entry!:
	unless(length($fasta_entry)<=2)
		{
		$fasta_entry =~ s/^>//;
		########################################
		process_entryF($fasta_entry);#
		########################################
		}



print "
seqs_kept:$seqs_kept
seqs not kept:$seqs_not_kept
";


my @seq_keys = keys  %store_sequences; @seq_keys = sort @seq_keys;

foreach my $id(@seq_keys)
	{
	print OUT1 ">$id\n$store_sequences{$id}";
	};

close OUT1;




};



##############################################################################################################




sub process_entryF
	{
	my $line = shift;
	my $current_id = "";
	if ($line =~ /^(.+)\n/ )
		{
		$current_id = $1;
		$line =~ s/^.+\n//;#$line =~ s/\n//;$line =~ s/\r//;my $current_seq_length = length($line);
		}else{
		die "\nerror 5631.\n"
		}


	if($current_id =~ /^[A-Z][a-z]+_[a-z]+$/)
		{
		$store_sequences{$current_id} = $line;
		$seqs_kept++;
		}else{
		print "removing non binomail:$current_id\n";
		$seqs_not_kept++;
		};



	}




##############################################################################################################













