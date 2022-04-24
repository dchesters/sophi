


# 
# 
# 
# 
# perl remove_fasta_entries.pl -fasta_in [file_name] -fasta_out [file_name] -rm_list [file_name] -regex [0/1] -remove/-retain -blastclust_list [0/1]
#
# 
#	format in remove_file, just a list:
#
# 	if($line =~ /(.+)/)
#		{my $rmID = $1;
# 
# change log
# 14 june 2014: in addition to looking for remove entries as full fasta id, 
#	option to look using a regex, useful if you want to remove all members of a species without specifying accessions
# 14 jan 2015: also retain fasta entries
# 05 Oct 2015: retain fasta entries works if user uses list of regex's
# 27 Jun 2019: gives error message if user inputs fasta file with Windows newlines
# 01 Nov 2019: options given on command line, 
# 12 Mar 2020: error message given if fasta file read where taxon list is expected
# 26 Nov 2020: prints file of representatives and all members
# 28 Nov 2020: prints list of anything couldnt find
# 01 Jul 2021: filters species of a comm table, command line option: -comm_table
# 22 Jul 2021: gives some screen output if run on huge slow file.
# 
# 
# 
# 
# 
# 
# 
# 
#############################################################################################################################



my $arg_string  = join ' ', @ARGV;


#####################################
read_command_arguments($arg_string);#
#####################################




# $database_file 		= $ARGV[0];
# $rm_file 		= $ARGV[1];
# $out_file 		= $ARGV[2];

# $regex 			= 0;	# default 0.
# $remove_or_retain	= 2;# 1=remove. 2=retain
# $blastclust_list 	= 1; # default = 0; 1= select one per OTU




##############################################################################################################



unless($database_file =~ /[\d\w]/ && $rm_file =~ /[\d\w]/ && $out_file =~ /[\d\w]/  )
	{
	die "\ncommand error, something missing.\n"
	};

if($blastclust_list == 1)
	{
	open (OTU_REPS, ">OTU_representatives_retained") || die "\nerror 66.\n";
	};


open(LIST_IN, $rm_file) || die "\nERROR 76. Cant open file:$rm_file.\n";
print "opened $rm_file\n";
while(my $line= <LIST_IN>)
	{
	$line =~ s/\n//;$line =~ s/\r//;
#	print "line:$line\n";

	if($line =~ /^>/){$looks_like_fasta++};

	if($blastclust_list == 1)
		{
		if($line =~ /^(\S+)\s+(\S.+)$/)
			{
			my $otu_representative = $1; my $other_members = $2;print OTU_REPS "$otu_representative\t$otu_representative $other_members\n";
			}elsif($line =~ /^(\S+)/)
			{
			my $otu_representative = $1; print OTU_REPS "$otu_representative\t$otu_representative\n";
			};
		$line =~ s/^(\S+)\s.*$/$1/;
		};

	if($line =~ /(.+)/)
		{
		my $rmID = $1; 
		 $lines_parsed++;
	#	if($lines_parsed <= 5)
	#		{
			print "$rmID ";
	#		};

		if($lines_parsed =~ /000$/){print "lines_parsed:$lines_parsed\trmID:$rmID\n"};
	#	 print "ID:($rmID)\n";

		#12S	Hylemya_variata_FJ025383
		#12S	Thambemyia_pagdeni_FJ808133

		$search_IDs{$rmID} = 1; # remove/retain list
		}
	}
close(LIST_IN);
close OTU_REPS;

print "
";

if($looks_like_fasta >= 20){print "\nWARNING, looks like fasta file has been read, where list is expected. maybe error\n"};


my @array_of_search_IDs = keys %search_IDs; # remove/retain list
@array_of_search_IDs = sort @array_of_search_IDs;


print scalar @array_of_search_IDs , " entries parsed, to be looked for in db for removal\n
";

#die;

######################################################################################



if($comm_table == 1)
	{
	print "\nfiltering a community ecology-type table. assuming taxonomic units are rows.\n";
	open(FASTA_IN, $database_file) || die "\nERROR 140. Cant open file $database_file.\n";
	open(OUT1, ">$out_file") || die"";

	my $linecount=0;	
	while(my $fasta_line= <FASTA_IN>)
		{
		$fasta_line =~ s/\n//;$fasta_line =~ s/\r//; # $fasta_line = $fasta_line . "\n";
		if($linecount == 0)
			{
			print OUT1 "$fasta_line\n";
			}else{
			my $print_current_entry =0;
			if($fasta_line =~ /^(\S+)/)
				{
				my $current_id = $1;

				if($search_IDs{$current_id} == 1)
					{
					#	print "$current_id is stored\n";
					if($remove_or_retain == 1) # REMOVE option
						{
						die "\nerror. option not implemented.\n";# $rm++; $lines_parsed2++;$print_current_entry=0;
						}else{ # RETAIN option
						$rt++;$print_current_entry=1; # print OUT1 ">$current_id\n$line";
						};
					};

				}else{
				print "warning, unexpected line of comm table:$fasta_line\n";
				};
			if($print_current_entry == 1){print OUT1 "$fasta_line\n";$comm_species_printed++};
			};
		$linecount++;
		};
	close FASTA_IN; close OUT1;
	print "comm_species_printed:$comm_species_printed\n";
	exit;	
	};






	open(FASTA_IN, $database_file) || die "\nERROR 184. Cant open file $database_file.\n";
	open(OUT1, ">$out_file") || die"";
	my $fasta_entry = "";
	while(my $fasta_line= <FASTA_IN>)
		{
		$fasta_line =~ s/\n//;$fasta_line =~ s/\r//;$fasta_line = $fasta_line . "\n";
		if($fasta_line =~ /^>.+/)
			{
			unless(length($fasta_entry)<=2)
				{
				$fasta_entry =~ s/^>//;
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

	unless(length($fasta_entry)<=2)
		{
		$fasta_entry =~ s/^>//;
		########################################
		process_entry($fasta_entry);#
		########################################
		}

	
print "finished printing to $out_file\n";
close OUT1;

open(COULDNT_FIND, ">ids_couldnt_find");
foreach my $id(@array_of_search_IDs)
	{
	if($fasta_ID_printed{$id} == 1)
		{	
		$ids_sucessfully_printed++;
		}else{
		$ids_couldnt_find++;  print COULDNT_FIND "$id\n";		
		};
	};
close COULDNT_FIND;

print "\n
 fasta_ID_count:$fasta_ID_count
 been removed from file:$rm
 retained:$rt
 ids_sucessfully_printed:$ids_sucessfully_printed
 ids_couldnt_find:$ids_couldnt_find [listed in file ids_couldnt_find]
\n\n";

unless($rt >= 1)
	{
	print "\nwarning, nothing retained. error can be caused by Windows newlines in fasta file.\n\n"
	};



print "

FIN.
";
exit;










###################################################################################################################

sub process_entry
{
my $line = shift;
my $current_id = "";
if ($line =~ /^(.+)\n/ )
	{
	$current_id = $1;
	$line =~ s/^.+\n//;#$line =~ s/\n//;$line =~ s/\r//;my $current_seq_length = length($line);
	}else{
	die "\nerror 5631.\n"
	};

$fasta_ID_count++;
if($fasta_ID_count < 5)
	{
	print "fasta ID:$current_id\n";
	}elsif($fasta_ID_count == 5)
	{
	print ".....\n\n"
	}elsif($fasta_ID_count =~ /0000$/)
	{
	print "$fasta_ID_count\n";	
	};

# print "current_id:$current_id\n";
	
my $print_current_entry;


if($regex == 1)
	{
	
	my $current_entry_matches_user_string=0;
	foreach my $search_regex(@array_of_search_IDs)
		{
		if($current_id =~ /$search_regex/)
			{
			#print "found entry to be removed ($current_id), matches remove regex ($search_regex)\n";
			$current_entry_matches_user_string =1;
			}
		}

	if($current_entry_matches_user_string == 1)
		{
		if($remove_or_retain == 1)# 1=remove. 2=retain
			{
			# current sequence matches string specified by user to be removed
			$rm++;$print_current_entry=0;
			}else{	
			# ... to be retained
			$rt++;$print_current_entry = 1; # print OUT1 ">$current_id\n$line"
			}
		
		}else{

		if($remove_or_retain == 1)
			{
			$rt++;$print_current_entry = 1; # print OUT1 ">$current_id\n$line";
			}else{
			$print_current_entry = 0;
			$rm++;
			};
		};


	}else{

# print "current_id2:$current_id\n";


 # my $current_idCOPY = $current_id;
 # $current_idCOPY =~ s/.+_(\d+)$/$1/;
 # if($rm_these{$current_idCOPY} == 1)	{$rm_these{$current_id}=1};

	if($search_IDs{$current_id} == 1)
		{
	#	print "$current_id is stored\n";
		if($remove_or_retain == 1) # REMOVE option
			{
			$rm++; $lines_parsed2++;$print_current_entry=0;
			}else{ # RETAIN option
			$rt++;$print_current_entry=1; # print OUT1 ">$current_id\n$line";
			}
		}else{
	#	print "$current_id not stored\n";
		if($remove_or_retain == 1)# REMOVE option
			{
			$rt++;$print_current_entry = 1; # print OUT1 ">$current_id\n$line";
			}else{
			$rm++;$print_current_entry = 0;
			};

		}
	# print "  print_current_entry:$print_current_entry\n";
	}



if($print_current_entry == 1){print OUT1 ">$current_id\n$line";$fasta_ID_printed{$current_id} = 1};



}

###################################################################################################################







#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




sub read_command_arguments
{

my $arguments = shift;
$arguments =~ s/\n//;$arguments =~ s/\r//;


if($arguments =~ /-fasta_in\s+(\S+)/)
	{
	$database_file = $1;
	}else{
	die "\ncommand error (-fasta_in)\n";
	};

if($arguments =~ /-fasta_out\s+(\S+)/)
	{
	$out_file = $1;
	}else{
	die "\ncommand error (-fasta_out)\n";
	};

if($arguments =~ /-rm_list\s+(\S+)/)
	{
	$rm_file = $1;
	}else{
	die "\ncommand error (-rm_list)\n";
	};

if($arguments =~ /-regex\s+(\d+)/)
	{
	$regex = $1;
	}else{
	die "\ncommand error (-regex)\n";
	};

if($arguments =~ /-remove/)
	{
	$remove_or_retain = 1;
	}elsif($arguments =~ /-retain/)
	{
	$remove_or_retain = 0;
	}else{
	die "\ncommand error (rm or retain)\n";
	};

if($arguments =~ /-blastclust_list\s+(\d+)/)
	{
	$blastclust_list = $1;
	}else{
	die "\ncommand error (-blastclust_list)\n";
	};

if($arguments =~ /-comm_table/)
	{
	$comm_table = 1;
	};

	
# $regex 			= 0;	# default 0.
# $remove_or_retain	= 2;# 1=remove. 2=retain
# $blastclust_list 	= 1; # default = 0; 1= select one per OTU






}#sub read_command_arguments



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################







