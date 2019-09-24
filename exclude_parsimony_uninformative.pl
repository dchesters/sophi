

#
#
#
#	exclude_parsimony_uninformative.pl 
#		perl script for processing fasta files alignment
#
#
#    	Copyright (C) 2015-2017  Douglas Chesters
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
#####################################################################################################
# 
# 
# 
# 
# need for this script identified by Ladislav Bocak
# for some applications , single taxon mutations can obscure pattern,
# only known function for removing variable but parsimony uninformative sites is in Paup,
# however, its severe limitation on number of characters excludes most modern day use.
# 
# these are the characters considered, 
# 	if $data_type == 1: [ACTG]
# 	elsif $data_type == 2) : [FLSYCWPHQRIMTNKSRVADEG]
# 	all other characters present are considered ambiguous
# thus, its important to correctly set dna or prot in the appropriate place below
# 
# 
# to run:
# perl exclude_parsimony_uninformative.pl input_file_name output_file_name
#
# requires fasta format input file. spaces (and other non [0-9a-zA-Z_-]) in fasta IDs are always a bad idea
# all characters except those listed above are ignored. if you need to consider other (e.g. YR) let me know
# if a column has two different states, and one of these only at one sequence,
# this column is uninformative and removed. 
# output written to file: parsimony_informative
# if you also want invarient or ambiguous sites removing, this is set in the options.
# 
# 
# tested on the following input: 
# >seq01
# AAAAAANAAAAAAAAA
# >seq02
# AAAANANAA?AAAAAA
# >seq03
# AAAAAAAAA?AATAAA
# >seq04
# AAAACACAAAAATT-N
# 
# result should be:
# >seq01
# AAAAAAAAAAAAA
# >seq02
# AAAAAAA?AAAAA
# >seq03
# AAAAAAA?AATAA
# >seq04
# AAAAAAAAAAT-N
# 
# 
# used for testing writing new partition files:
# >tax1
# AAAAAA??AAAAAAAAAAAATTTTNNTTTTCCCCCGG
# >tax2
# AAAAAAA???AAAAAAAAAATTTTNTTTTTCCCCCGG
# 
# DNA, part1 = 1-20
# DNA, part2 = 21-30
# DNA, part3 = 31-35
# DNA, part4 = 36-37
# 
# 
# 
# 
# ###################################################################################################
# 
# 
# 
# 	CHANGES
# 	2015-11-11: quicker implementation, previously was too slow on genomic matrices.
# 			also, option for removing invariant sites also (again ignoring ambiguous states)
# 			also, can use on protein seqs, just set datatype.
# 	2015-11-12: some minor issues corrected
# 			option to remove columns of only ambiguous states 
# 			these are not insignificant when using a typical large processed phyloinformatic matrix
# 	2017-11-05: takes partition file as input, and makes new one adjusted for removed columns
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
#############################################################################################################





$arg_string  = join ' ', @ARGV;


#####################################
read_command_arguments($arg_string);#
#####################################


# $input 		= $ARGV[0];
# $output 	= $ARGV[1];


#############################################################################################################


	#      OPTIONS:

# $data_type 			= 2; #1=DNA; 2=PROT
# $remove_invariant_columns 	= 1; # while we are at it, remove other types of useless columns ...
# $remove_ambiguous_columns 	= 1; # if =1, remove columns composed only ?- etc .. 
						# these can occur due to post-alignment processing, 
					 	



#############################################################################################################




unless($input =~ /[\d\w]/ && $output =~ /[\d\w]/)
	{die "\ncommand error.\nto run:\nperl exclude_parsimony_uninformative.pl input_file_name output_file_name\n\n"};
if($input eq $output )
	{die "\ncommand error: input eq output\n\n"};

print "\nexclude_parsimony_uninformative.pl
input:$input
output:$output\n";

open(IN_FILTER , "$input") || die "\nerror. cant open\n";
open(OUT, ">$output") || die "\nerror cant open outfile\n";
open(OUT2, ">parsimony_informativeLOG") || die "\nerror cant open outfile\n";




#############################################################################################################



	# READ FILE:

my $file_as_string = "";
print "reading file $input\n";
while (my $line = <IN_FILTER>)	
	{
	$file_as_string .= $line;my $file_length = length($file_as_string);
	if($line =~ /^>(.+)/)
		{
		my $id = $1;$fasta_count++;
		#print "fasta entry:$id, total count:$fasta_count, file size:$file_length\n";
		};
	};close(IN_FILTER);
print "\nfile is read, $fasta_count entries\n";
my @all_lines = split />/, $file_as_string;
my $count_entries  = $#all_lines;
unless($#all_lines == $fasta_count){die "\nsome kind of error trying to read fasta entries, check format\n"}
print OUT2 $#all_lines , " seqs in file\n";



#############################################################################################################



	# get sequence length:

my $alignment_length= "";
my $line = $all_lines[1];
if($line =~ /^(.+)/)
	{
	$line =~ s/^.+\n//;$line =~ s/\012\015?|\015\012?//g;$line =~ s/\n//g;$line =~ s/\r//g;$line =~ s/[\s\t]//g;
	$alignment_length = length($line);
	}else{die "\nerror 23\n"};
print OUT2 "alignment_length:$alignment_length\n";
print "alignment_length:$alignment_length\n";



#############################################################################################################


# for genomic data its too slow to , for each column read through whole alignment just to get that column,
# so need to put each column into some instant access form first:

print "storing columns...\n";
#for my $j(0 .. ($alignment_length - 1))
#	{
	#if($j =~ /00$/){print "column $j of $alignment_length\n"};
	#my @current_column = ();
	for my $each_line(1 .. $#all_lines)# each fasta entry
		{
		my $line = $all_lines[$each_line];#print "\nline:$line\n";
		if($line =~ /^(.+)/)
			{
			my $speciesid = $1;	#print "$speciesid\n";
			$line =~ s/^.+\n//;$line =~ s/\012\015?|\015\012?//g;$line =~ s/\n//g;$line =~ s/\r//g;$line =~ s/[\s\t]//g;
			}
		for my $j(0 .. ($alignment_length - 1))
			{
			my $char = substr $line , $j , 1;
			$columns{$j} .= $char;
			};
		#push @current_column , $char;

		if($#all_lines >= 10000 && $each_line =~ /0000$/){print "$each_line of $#all_lines\n"};

		}

#	};
print "columns stored for quick access\n";



#############################################################################################################



	# look through each column:

print "look through each column ...\n";

my @record_parsimony_uniformative_columns = ();

for my $j(0 .. ($alignment_length - 1))
	{
	if($j =~ /00000$/){print "column $j of $alignment_length\n"};
	my @current_column = ();
	@current_column = split // , $columns{$j};
	print OUT2 "\ncolumn " , ($j+1) , "\n\tchars:@current_column\n";
	#print "\ncolumn " , ($j+1) , "\n\tchars:@current_column\n";

	#my $count_nucleotides =0;
	my %current_column_bases  =();

	for my $base(@current_column)
		{
	#	print "\tbase:$base\n";
		if($data_type == 1)
			{
			if($base =~ /[ACTG]/i)
				{
				$current_column_bases{$base}++;
				$char_present++;
				}elsif($base =~ /[nN\?\-]/)
				{
				$char_absent++;
				}else{
				$char_not_read++;
				};
			}elsif($data_type == 2)
			{
			# these are from NCBI's 'standard' code, others are available, e.g. invert mt
			if($base =~ /[FLSYCWPHQRIMTNKSRVADEG]/i)
				{
				$current_column_bases{$base}++;
				$char_present++;
				}elsif($base =~ /[\?\-]/)
				{
				$char_absent++;
				}else{
				$char_not_read++;
				};
			}else{
			die "\nerror setting data type\n";
			}
		}



	my @current_bases_keys  =keys %current_column_bases;
	my $number_of_bases_current_column = scalar @current_bases_keys;
	my $current_column_has_state_for_just_one_sequence_only=0;
	foreach my $base(@current_bases_keys)
		{
		my $count = $current_column_bases{$base};
		print OUT2 "\t\tbase:$base count:$count\n";
		# correct me if wrong, as far as i can see, 
		# parsimony uniformative will be those with two different bases in the column,
		# and one of those bases will be present only in one sequence:
		if($count == 1){$current_column_has_state_for_just_one_sequence_only = 1};
		}

	my $currnet_column_is_parsimony_uniformative  =0;
	if($current_column_has_state_for_just_one_sequence_only == 1 && 
		$number_of_bases_current_column == 2)
			{
			print OUT2 "\tPARSIMONY UNINFORMATIVE\n";
			$currnet_column_is_parsimony_uniformative = 1;
			$count_parsimony_uninformative_columns++;
			$count_rm_columns++;
			};

	if($number_of_bases_current_column == 1)
		{
		$count_invariant_columns++;
		print OUT2 "\tINVARIANT\n";
		if($remove_invariant_columns == 1)
			{$currnet_column_is_parsimony_uniformative = 1;
			$count_rm_columns++}
		}

	if($number_of_bases_current_column == 0)
		{
		$count_ambiguous_columns++;
		print OUT2 "\tAMBIGUOUS\n";
		if($remove_ambiguous_columns == 1)
			{$currnet_column_is_parsimony_uniformative = 1;
			$count_rm_columns++}
		}

print OUT2"\trm col:$currnet_column_is_parsimony_uniformative total rmd cols:$count_rm_columns
	total PI:$count_parsimony_uninformative_columns
	total INV:$count_invariant_columns
	total AMB:$count_ambiguous_columns
";
	$record_parsimony_uniformative_columns[$j] = $currnet_column_is_parsimony_uniformative;

	$count_columns++;
	}; # for my $j(0 .. ($alignment_length - 1))

###########################################################################################################

$total_read = $char_present + $char_absent + $char_not_read;
$data_presence = $char_present / $total_read;

print "
char_present:$char_present
char_absent:$char_absent
char_not_read:$char_not_read
data_presence:$data_presence
";



	# finally, go through seqs again, print all except the pars uninfom sites:

for my $each_line(1 .. $#all_lines)
	{
	my $line = $all_lines[$each_line];#print "\nline:$line\n";
	my $speciesid;
	if($line =~ /^(.+)/)
		{
		$speciesid = $1;	#print "$speciesid\n";
		$line =~ s/^.+\n//;$line =~ s/\012\015?|\015\012?//g;$line =~ s/\n//g;$line =~ s/\r//g;$line =~ s/[\s\t]//g;
		}
	print OUT ">$speciesid\n";

	for my $j(0 .. ($alignment_length - 1))
		{
		my $char = substr $line , $j , 1;
		if($record_parsimony_uniformative_columns[$j] == 0)
			{
			print OUT "$char";
			}
		};
	print OUT "\n";
	}

close OUT;




##############################################################################################################


##############################################################################################################



# if re-printing partition file

if($reprint_partition_file == 1)
{
print "\nprinting new partition files with site numbers recalculated\n";

$previous_partition_name;

$new_alignment_position = 0;

for my $j(0 .. ($alignment_length - 1))
	{
	my $is_column_PI = $record_parsimony_uniformative_columns[$j];
	# site_is_from_partition_name starts at index 1	, whereas record_parsimony_uniformative_columns is index 0
	my $adjust_index = $j+1;my $part_name = $site_is_from_partition_name{$adjust_index};
	unless($part_name =~ /[\d\w]/){die "\nerror , no partition name found for site $j / $adjust_index\n"};

#	print "j:$j is_column_PI:$is_column_PI part_name:$part_name\n";


	if($is_column_PI == 0)
		{
		if(exists($store_partition_names{$part_name}))
			{
		
			}else{
		#	print "new partiiton encountered, named $part_name\n";
			$new_start_position{$part_name} = $new_alignment_position+1;
			$new_end_position{$previous_partition_name} = $new_alignment_position;
			$previous_partition_name = $part_name;
			push @new_partitions, $part_name;
			};
		$store_partition_names{$part_name} = 1;
		$new_alignment_position++;
		}; # if($is_column_PI == 0)


#	if($j >= 20){die ""};
	};
$new_end_position{$previous_partition_name} = $new_alignment_position;

open(NEW_PARTITIONS1, ">new_partitions.rax") || die "\nerror 359\n";
open(NEW_PARTITIONS2, ">new_partitions.cs") || die "\nerror 359\n";
open(NEW_PARTITIONS3, ">new_partitions.pf") || die "\nerror 359\n";

foreach my $partition(@new_partitions)
	{
	my $start  = $new_start_position{$partition};my $end = $new_end_position{$partition};

	# raxml format
	# PROT, gene4_insectaNUCL_EOG505QG6_clo_pruned_fa_RMmulti = 4374-5926
	print NEW_PARTITIONS1 "$dtype, $partition = $start-$end\n";

	# charset
	# charset gene5_insectaNUCL_EOG508KQM_clo_pruned_fa_RMmulti = 5927 - 6674 ;
	print NEW_PARTITIONS2 "charset $partition = $start - $end ;\n";

	# partiiton finder
	# gene5_insectaNUCL_EOG508KQM_clo_pruned_fa_RMmulti = 5927-6674;
	print NEW_PARTITIONS3 "$partition = $start-$end;\n";

	};
close NEW_PARTITIONS1;
close NEW_PARTITIONS2;
close NEW_PARTITIONS3;

};# if($reprint_partition_file == 1)

##############################################################################################################


##############################################################################################################







# dont print this if using genomic alignment, the log file probably wont be openable
if($#record_parsimony_uniformative_columns <= 10000)
	{
	print OUT2 "\nfollowing columns to be removed:\n@record_parsimony_uniformative_columns\n";
	};

$new_alignment_length = $count_columns - $count_rm_columns;

print "
count_columns:$count_columns
count_parsimony_uninformative_columns:$count_parsimony_uninformative_columns
count_invariant_columns:$count_invariant_columns
count_ambiguous_columns:$count_ambiguous_columns
count_rm_columns:$count_rm_columns
new_alignment_length:$new_alignment_length
";


print "
	exclude_parsimony_uninformative.pl
	processed alignment written to $output
	further details written to parsimony_informativeLOG

FIN.

";






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

print "";


if($arguments =~ /-in\s+(\S+)/)
	{
	$input = $1;
	open(INTEST , $input) || die "\ninput not found in current directlry:$input\n";
	close INTEST;
	}else{
	print "\nerror reading command arguments, -in input\n";die""
	};
if($arguments =~ /-out\s+(\S+)/)
	{
	$output = $1;
	}else{
	print "\nerror reading command arguments, -out output\n";die""
	};
if($arguments =~ /-data_type\s+prot/i)
	{
	$data_type = 2;
	}elsif($arguments =~ /-data_type\s+dna/i)
	{
	$data_type = 1;
	}else{
	print "\nerror reading command arguments, -data_type [prot/dna]\n";die""
	};
if($arguments =~ /-remove_invar/)
	{
	$remove_invariant_columns 	= 1;
	}else{
	$remove_invariant_columns 	= 0;
	};
if($arguments =~ /-remove_ambig/)
	{
	$remove_ambiguous_columns 	= 1;
	}else{
	$remove_ambiguous_columns 	= 0;
	};

if($arguments =~ /-partition_file\s+(\S+)/)
	{
	my $partition_filename = $1;$reprint_partition_file = 1;
	open(PART, $partition_filename) || die "\nerror. user specified partition file named $partition_filename, but cannot open it.\n";
	while (my $line = <PART>)
		{
		$line =~ s/\n//;$line =~ s/\r//;
		# PROT, gene3_insectaNUCL_EOG502V81_clo_pruned_fa_RMmulti = 3802-4373
		if($line =~ /^(\S+)\,\s(\S+)\s\=\s(\d+)\-(\d+)$/)
			{
			$dtype = $1; my $part_name = $2; my $part_start = $3; my $part_fin = $4;
		#	print "\nline:$line\n\tname:$part_name start:$part_start end:$part_fin\n";

			for my $site_number($part_start .. $part_fin)
				{
				$site_is_from_partition_name{$site_number} = $part_name; # print "$site_number ";
				};

			$count_partitions_parsed++;			
			}elsif($line =~ /[\w\d]/)
			{
			print "warning, cant parse this line of partition file:$line\n";
			};
		};
	close PART;
	print "\nuser specified partition file named $partition_filename.\n\tfound $count_partitions_parsed\n\n";

	unless($count_partitions_parsed >= 2){die "\nerror, too few partitions read from file \n"}
	}else{

	};


# perl ~/usr_scripts/alignment_processing/exclude_parsimony_uninformative.pl 
# -in insectaNUCL.smatrix.RM2 -out insectaNUCL.smatrix.RM2.epu -data_type prot -remove_invar -remove_ambiguous

# $input 		= $ARGV[0];
# $output 	= $ARGV[1];
# $data_type 			= 2; #1=DNA; 2=PROT
# $remove_invariant_columns 	= 1; # while we are at it, remove other types of useless columns ...
# $remove_ambiguous_columns 	= 1; # if =1, remove columns composed only ?- etc .. 





};






