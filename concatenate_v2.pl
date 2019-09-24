
###################################################################################################################################
#
#
#
#	concatenate_v2.pl 
#	perl script for concatenating dna sequence data from different fasta files
#
#
#    	Copyright (C) 2013-2018  Douglas Chesters
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
###################################################################################################################################
#
#
#
#	to run:
#	perl concatenate_v2.pl -missing_data_char ? -remove_accession 2 -matrices sequences1.fas sequences2.fas ...........
#
#
#	change log
#	01 Oct 2013: added gpl
#	05 Aug 2014: bugfix on printing partition file
#	30 Sep 2014: warnings given regarding accession removal
#	27 Jan 2015 (v2.03): removed accession only from binomials
#	01 Apr 2015: option to print only concatenates which have sequences of user specified genes.
#	09 Nov 2015: outputs two partition description files, 
#		one in raxml's format, one 'charset style' which other software uses (e.g. MARE)
#		charset YourFirstPartition = 1 - 100 ;
#	12 Nov 2015: some minor issues when concatenating protein sequences.
#			and trims the accession delimter used by phylotreepruner '|'
#	09 Apr 2016: all options are given by user in command instead of written in script
#			more suitable as part of a pipeline
#	21 May 2016: minor bugfix in command argumnents
#	19 May 2017: option to trim accession delimited by '|' character (phylotreepruner output)
#	20 Oct 2017: prints error message when incorrect number of partitions is specified in -required_data option
#	05 Nov 2017: prints partition finder format datablock
#		bugfix first partition was not defined in one of the output files
#		option to force data type in partition file; argument -data_type with options DNA, PROT or BIN
#		data type autodetect improved marginally, although user definition is preferable
#	11 Jan 2018: in case you are doing phylogenomics, have already done meta-partitioning  (by partition_finder)
#		but want to use these on a set of orthologs subsequently modified ....
#		takes the previous set of meta-partitions,
#		then writes the overlapping set out for the regular single gene input files
#	31 Mar 2018: prints first gene (usually coi) with ids as used in supermatrix
#	
#	
#
#
#
#
#
#
#
######################################################################################################



my $arg_string  = join ' ', @ARGV;

#####################################
read_command_arguments($arg_string);#
#####################################



# @input_array = @ARGV;unless(scalar @input_array >= 2){die "\n\nError, did you specify input alignment? (concatenate v2)\n\n"};


			# warning, if you are accustomed to DNA seqs and start using proteins, 
			# you cannot use N as ambiguous char, it is already used as a character state!
# $missing_data_char = "?";######################################################

# $remove_accession 	= 2;	# assumed format for fasta sequence identifiers: species id, followed by underscore, followed by accession
				# 0 = dont RM
				# 1 = remove terminal string, s/_[^_]+$//;
				# 2 = strict setting, remove terminal string only from binomials
				#	s/^([A-Z][a-z]+_[a-z]+)_[^_]+$/$1/;

# $genus_level		= 0; # default 0. if 1, will trim of species namse if present, and concat at the genus level

# $required_data = ""; # default, nothing in this object. 
					# string like this "Y--", means you are concatenating 3 genes
					# but you only want to print a concatenate if it possesses a seq of the first gene
					# NOTE: use Y for reqruired gene for making a concatenate, 
					# otherwise use your missing data character for gene not neccessary
					# i.e. not literally '-'. 


@datatypes =();
@gene_lengths = ();




my %species_printed;


################################
concatenate_these_alignments();#
################################


########################
print_partition_file();#
########################


if ($#input_array <= 24)
	{
	##############################
	print_processed_alignments();#
	##############################
	};


if($remove_accession == 1)
	{
	print "\nWarning. \$remove_accession is set to 1. 
this is undiscriminatory , and only to be used where strict binomials are used in all cases
if in doubt, set to 0 or 2
\n";
	};
if($remove_accession == 2)
	{
	print "
remove accesion option set to 2,
this will work on either 
Genus_species_accession
Genus_species|accession
latter used by phylotreepruner
note there might be problems if you have both _ and | in your IDs ...

"	};
if($remove_accession == 3)
	{
	print "
remove accesion option set to 3,
this will work only on 
Genus_species|accession
 used by phylotreepruner
";
	};

print "end of script\n";
exit;






########################################################################################################################


sub concatenate_these_alignments
{

my %sequence_hash = ();
my %species_hash = ();

print "\n\nsub concatenate_these_alignments\n";


my $locus_count = 0;
my @missing_data_array = ();
my %duplicate_id_count_hash=();


foreach my $current_processed_alignment(@input_array)
	{
	print "alignment:$current_processed_alignment\n";


	open(CURRENT_PROCESSED_ALIGNMENT , $current_processed_alignment) || die "\n\nerror cant open $current_processed_alignment\n\n";	
	my $file_as_string9 = "";
	while (my $line77 = <CURRENT_PROCESSED_ALIGNMENT>)
		{$file_as_string9 .= $line77;#print 
		}
	close(CURRENT_PROCESSED_ALIGNMENT);
	my @all_lines9 = split />/, $file_as_string9;


	print scalar @all_lines9 , " seqs in file\n";


# need to do this after the species to be used, have been selected
#open(CURRENT_ALIGNMENT , ">$current_processed_alignment.accession_rm") || die "\n\nerror cant open $current_processed_alignment\n\n";	


	for my $each_line7(1 .. $#all_lines9)
		{
		my $line7 = $all_lines9[$each_line7];

		if($line7 =~ /(.+)\n/)
			{
			my $speciesid = $1;
			if($speciesid =~ /^[A-Z][a-z]+_[a-z]+_[^_]+$/ && $remove_accession == 0)
				{print "\nWARNING. you have set remove accession as 0, " . 
				"but it looks like there are accessions in your lables. concatenation will not work properly.\n"}
			if($remove_accession == 1)
				{
				$speciesid =~ s/_[^_]+$//;
				};
			if($remove_accession == 2)
				{
				$speciesid =~ s/^([A-Z][a-z]+_[a-z]+)[_\|][^_\|]+$/$1/;
				};
			if($remove_accession == 3)
				{
				$speciesid =~ s/\|[^\|]+$//;
				};

			if($genus_level == 1)
				{
				#$speciesid =~ s/^([^7]+7).+/$1/; # no longer using those tobycodes
				$speciesid =~ s/^([A-Z][a-z]+)_.+$/$1/; 
				};


		#	$speciesid =~ s/BOLD.+/$1/;

			$line7 =~ s/^.+\n//;
			$line7 =~ s/\n//g;$line7 =~ s/\r//g;$line7 =~ s/\012\015?|\015\012?//g;

			#print CURRENT_ALIGNMENT ">$speciesid\n$line7\n";

			$sequences{$speciesid . "______" . $current_processed_alignment} = $line7;#print "line7:$line7";die;


			# this is just used when writing partition files:

			if($force_datatype =~ /\w/)
				{
				$datatypes[$locus_count] = $force_datatype;
				}else{
				

				if($line7 =~ /[actgACGT]{9}/)
					{
					$datatypes[$locus_count] = "DNA"
					}elsif($line7 =~ /[FLSYCWPHQRIMTNKSRVADEG]{3}/i)
					{
					$datatypes[$locus_count] = "PROT"
					}else{
					$datatypes[$locus_count] = "BIN"
					};

				}


			my $hashkey = $speciesid . "______" . $current_processed_alignment;

			$duplicate_id_count_hash{$hashkey}++;

			if($duplicate_id_count_hash{$hashkey}>=2)
				{
				print "$hashkey , $duplicate_id_count_hash{$hashkey} >= 2\n"
				}else{
			# problem writing first seqeunce observed for species, instead of last which is used in concat
			#	if($locus_count == 0){print FIRST_GENE ">$speciesid\n$line7\n"};
				};

			$species_hash{$speciesid} = 1;

			# sequence is stored already, this just makes an object containing 
			# apporiate size of missing data for each gene
			$line7 =~ s/\S/$missing_data_char/g;
			$missing_data_array[$locus_count] = $line7;$line7 =~ s/\s//g;$gene_lengths[$locus_count] = length($line7);

			}
		}
	$locus_count++;
	

#close CURRENT_ALIGNMENT;

	}


open(FIRST_GENE, ">first_gene_processed_IDs") || die "";


my @allsp = keys %species_hash;
@allsp = sort @allsp;
$number_entries = scalar @allsp;

print "\n\nalignments read. now printing supermatrix\n";

open(OUT12 , ">current_supermatrix2") || die "\nerror 1168\n";

my $count_concatenates = 0;
my %duplicate_ids=();
my %duplicate_ids2=();
my %count_partition_patterns;

foreach my $sp(@allsp)
	{

	#print OUT12 ">$sp\n";
	my $print_string = ">$sp\n";

	$locus_count = 0;
	print "$sp\n";

	my $partition_pattern = "";
	foreach my $current_processed_alignment(@input_array)
		{
		my $hash_key = $sp . "______" . $current_processed_alignment;

		if(exists($sequences{$hash_key}))
			{
			#print OUT12 $sequences{$hash_key};
			$print_string .= $sequences{$hash_key};
			if($locus_count == 0){print FIRST_GENE ">$sp\n$sequences{$hash_key}\n"};

			if($duplicate_id_count_hash{$hash_key}>=2)
				{
				$duplicate_ids{$current_processed_alignment}++;
				$duplicate_ids2{$current_processed_alignment}+= ($duplicate_id_count_hash{$hash_key}-1);
				
				}
			$partition_pattern .= "Y";
			}else{
		#	print "not found\n";
			#print OUT12 $missing_data_array[$locus_count];
			$print_string .= $missing_data_array[$locus_count];

			$partition_pattern .= $missing_data_char; #"N"; 2015-11-11: no good for AAs, oops
			}			# although this only comes into spaly if requried partition patterns are set
		$locus_count++;
		}	
	$count_partition_patterns{$partition_pattern}++;
	#print OUT12 "\n";	
	$print_string .= "\n";


	my $print_current_line = 1;
	if($required_data =~ /./)
		{
		my @split1 = split // , $required_data;
		my @split2 = split // , $partition_pattern;
		unless($#split1 == $#split2){die "\nerror 251. check script settings \$remove_accession and \$required_data \n"};

		foreach my $index(0 .. $#split1)
			{	
			my $test1 = $split1[$index];
			my $test2 = $split2[$index];
			#print "test1:$test1 test2:$test2\n";
			if($test1 eq "Y" && $test2 eq $missing_data_char){$print_current_line = 0}
			};

		}else{
		$print_current_line = 1;
		};

	if($print_current_line == 1)
		{
		print OUT12 $print_string;
		$species_printed{$sp}=1;
		$count_concatenates++;
		}


	}#foreach my $sp(@allsp)

close(OUT12);
close FIRST_GENE;



print "\ncount_concatenates:$count_concatenates\n\n";

my @keys5 = keys %duplicate_ids;@keys5 = sort @keys5;

for my $k(@keys5)
	{
#	print "number of species IDs that have more than 1 motus:$duplicate_ids{$k} in $k\n";
#	print "number of motus that are in the duplicated class: $duplicate_ids2{$k} in $k\n\n";
	}

open(OUTINFO , ">partition_patterns") || die "\nerror 253\n";

foreach my $current_processed_alignment(@input_array)
	{
	print OUTINFO "$current_processed_alignment\t";
	};
print OUTINFO "\n";

my @partit_patterns = keys %count_partition_patterns;@partit_patterns = sort @partit_patterns;
foreach my $pattern(@partit_patterns)
	{
	my $count = $count_partition_patterns{$pattern};
	print OUTINFO "$pattern\t$count\n";
	};

close OUTINFO;

};#sub concatenate_these_alignments






##########################################################################################################################
#
#
#
#
##########################################################################################################################








sub print_partition_file
{


open(OUT_PF , ">partitionfile2") || die "cant open file\n";
open(OUT_CF , ">charsetfile2") || die "cant open file\n";

open(OUT_PARTITFINDER , ">partitionfinder_datablock") || die "cant open file\n";


$current_start =1;
$current_end = $gene_lengths[0];
my $curfile_name = "gene0" . "_" . $input_array[0];$curfile_name =~ s/[^\d\w]/_/g;
print OUT_PF $datatypes[0], ", " , $curfile_name , " = " , $current_start ,"-", $current_end , "\n";
print OUT_CF "charset $curfile_name = $current_start - $current_end ;\n";
print OUT_PARTITFINDER $curfile_name , " = " , $current_start ,"-", $current_end , ";\n";
	my $modified_id = $input_array[0];$modified_id =~ s/(RMmulti)_order/$1/;
	$modified_id =~ s/[^\d\w]/_/g;
	my $start_and_end = $current_start . "-" . $current_end;$start_and_ends{$modified_id} = $start_and_end;



# partit_belongs_to_metapartit
# model_of_metapartition

# charset gene1_insectaNUCL_EOG50001C_clo_pruned_fa_RMmulti_order = 1746 - 2936 ;
#   PROT, gene1_insectaNUCL_EOG50001C_clo_pruned_fa_RMmulti_order = 1746-2936

# metapartition style:
# gene565_insectaNUCL_EOG5W6MBD_clo_pruned_fa_RMmulti, gene153_insectaNUCL_EOG579CQD_clo_pruned_fa_RMmulti, gene333_insectaNUCL_EOG5H9W25_clo_pruned_fa_RMmulti, gene2_insectaNUCL_EOG502V7Z_clo_pruned_fa_RMmulti



print "\tcurrent_start:$current_start current_end:$current_end\n";

for $i(1 .. $#input_array)	# for each file name
	{

	$current_start = $current_end+1;
	$current_end = $current_end+$gene_lengths[$i];
	my $curfile_name = "gene$i" . "_" . $input_array[$i];$curfile_name =~ s/[^\d\w]/_/g;

	print OUT_PF $datatypes[$i], ", " , $curfile_name , " = " , $current_start ,"-", $current_end , "\n";
	print OUT_PARTITFINDER $curfile_name , " = " , $current_start ,"-", $current_end , ";\n";
	#COI_2ndpos = 651-1415\3;
	#COI_3rdpos = 652-1415\3;
	#16S = 1416-1897;

	#charset YourFirstPartition = 1 - 100 ;
	print OUT_CF "charset $curfile_name = $current_start - $current_end ;\n";

	my $modified_id = $input_array[$i];$modified_id =~ s/(RMmulti)_order/$1/;
	$modified_id =~ s/[^\d\w]/_/g;
	my $start_and_end = $current_start . "-" . $current_end;
	$start_and_ends{$modified_id} = $start_and_end;
	print "stpring id:$modified_id\n";

#	print "i:$i gene_lengths[i]:$gene_lengths[$i]\n\tstart:$current_start end:$current_end\n";

# raxml style:
#JTTF, Subset326 = 250378-250628, 293559-293919, 286582-286954, 307328-307778, 274494-274763


    	}

close(OUT_PF);
close(OUT_CF);
close (OUT_PARTITFINDER);

print "
" , scalar @input_array , " partitions printed\nnumber_fasta entries:$number_entries\n";




#######################################################################################################



	# if re-writing meta-partitions


if($print_meta_parts == 1)
{
print "\nprinting adjusted meta partitions\n";
my @list_mps = keys %meta_partitions;@list_mps = sort @list_mps;

open(NEW_METAPARTS , ">new_metapartitions");
close NEW_METAPARTS;
open(NEW_METAPARTS_INFO , ">new_metapartitions_INFO");
close NEW_METAPARTS_INFO;


foreach my $metapart_number(@list_mps)
	{
	my $meta_parttiion = $meta_partitions{$metapart_number};# these are the ORIGINAL meta-partitions
									# i.e. will sometimes have orthologs which arnt in the new data
	my @split_partition = split /\,/ , $meta_parttiion;
	my @current_new_metaparttion = ();my @current_new_metaparttion_positions = ();
	my $model = $model_of_metapartition{$metapart_number};
	my $raxml_model_code = $raxml_codes{$metapart_number};

	foreach my $part(@split_partition)
		{
		# :gene414_insectaNUCL_EOG5N5TCD_clo_pruned_fa_RMmulti
		$part =~ s/^gene\d+_(insectaNUCL_.+_clo_pruned_fa_RMmulti)$/$1/;
		# $partit_belongs_to_metapartit{$part} = $metapart_number
		print "metapart_number:$metapart_number part:$part\n";

		if(exists($start_and_ends{$part}))
			{
			push @current_new_metaparttion_positions , $start_and_ends{$part};
			push @current_new_metaparttion , $part;
			}else{
		#	die "\nerror 5988999\n";
			};
		};

	if(scalar @current_new_metaparttion >= 1)
		{
		my $print_this = join ', ' , @current_new_metaparttion_positions;
		my $print_this_geneIDs = join ', ' , @current_new_metaparttion;
		print "\t$print_this\n";

		open(NEW_METAPARTS , ">>new_metapartitions");
		print NEW_METAPARTS "$raxml_model_code, metapart$metapart_number = $print_this\n";
		close NEW_METAPARTS;

		open(NEW_METAPARTS_INFO , ">>new_metapartitions_INFO");
		print NEW_METAPARTS_INFO "$metapart_number $model $raxml_model_code $print_this_geneIDs $print_this\n";
		close NEW_METAPARTS_INFO;


		};

	
	};



};


#######################################################################################################





};# sub print_partition_file



##########################################################################################################################
#
#
#
#
##########################################################################################################################




sub print_processed_alignments
{



foreach my $current_processed_alignment(@input_array)
	{
	print "alignment:$current_processed_alignment\n";
	my $species_ARE_printed_to_processed_single_gene_alignment =0;
	my $species_NOT_printed_to_processed_single_gene_alignment =0;


	open(CURRENT_PROCESSED_ALIGNMENT , $current_processed_alignment) || die "\n\nerror cant open $current_processed_alignment\n\n";	
	my $file_as_string9 = "";
	while (my $line77 = <CURRENT_PROCESSED_ALIGNMENT>)
		{$file_as_string9 .= $line77;#print 
		}
	close(CURRENT_PROCESSED_ALIGNMENT);
	my @all_lines9 = split />/, $file_as_string9;
	print "\t" , scalar @all_lines9 , " seqs in file\n";


	#open(CURRENT_ALIGNMENT , ">$current_processed_alignment.accession_rm") || die "\n\nerror cant open $current_processed_alignment\n\n";	


	for my $each_line7(1 .. $#all_lines9)
		{
		my $line7 = $all_lines9[$each_line7];

		if($line7 =~ /(.+)\n/)
			{
			my $speciesid = $1;
			if($speciesid =~ /^[A-Z][a-z]+_[a-z]+_[^_]+$/ && $remove_accession == 0)
				{print "\nWARNING. you have set remove accession as 0, " . 
				"but it looks like there are accessions in your lables. concatenation will not work properly.\n"}
			if($remove_accession == 1)
				{
				$speciesid =~ s/_[^_]+$//;
				};
			if($remove_accession == 2)
				{
				$speciesid =~ s/^([A-Z][a-z]+_[a-z]+)_[^_]+$/$1/;
				};

			$line7 =~ s/^.+\n//;
			$line7 =~ s/\n//g;$line7 =~ s/\r//g;$line7 =~ s/\012\015?|\015\012?//g;

			if(exists($species_printed{$speciesid}))
				{
				print CURRENT_ALIGNMENT ">$speciesid\n$line7\n";
				$species_ARE_printed_to_processed_single_gene_alignment++;
				}else{
				$species_NOT_printed_to_processed_single_gene_alignment++;
				
				}
			}
		};

print "\tspecies_ARE_printed_to_processed_single_gene_alignment:$species_ARE_printed_to_processed_single_gene_alignment
\tspecies_NOT_printed_to_processed_single_gene_alignment:$species_NOT_printed_to_processed_single_gene_alignment
";

	};#foreach my $current_processed_alignment(@input_array)





}

##########################################################################################################################





sub read_command_arguments
{

my $arguments = shift;
$arguments =~ s/\n//;$arguments =~ s/\r//;


# @input_array = @ARGV;

if($arguments =~ /\-matrices\s+(.+)/)
	{
	my $list = $1;
	@input_array = split /\s+/ , $list;
	unless(scalar @input_array >= 2){die "\n\nError, did you specify input alignment? (concatenate v2)\n\n"};
	}else{
	die "\ncommand error, this is missing:\n-matrices matrix1 matrix2 matrix3";
	}

			# warning, if you are accustomed to DNA seqs and start using proteins, 
			# you cannot use N as ambiguous char, it is already used as a character state!
# $missing_data_char = "?";######################################################
if($arguments =~ /\-missing_data_char\s+(\S+)/)
	{	
	$missing_data_char = $1;
	}else{
	$missing_data_char = "?";
	};

if($arguments =~ /\-remove_accession\s+(\S+)/)
	{	
	$remove_accession = $1;
	}else{
	$remove_accession = 0;
	};
# $remove_accession 	= 2;	# assumed format for fasta sequence identifiers: species id, followed by underscore, followed by accession
				# 0 = dont RM
				# 1 = remove terminal string, s/_[^_]+$//;
				# 2 = strict setting, remove terminal string only from binomials
				#	s/^([A-Z][a-z]+_[a-z]+)_[^_]+$/$1/;
				# 3 = trim after |, which will be phylotreepruner output


if($arguments =~ /\-metapartitions\s+(\S+)/)
	{	
	$metapartitions = $1;$print_meta_parts = 1;

	print "
you have specified to read and re-write metapartitions, 
relevent partition finder file is usually called best_scheme.txt
your input:$metapartitions,
";

	open(IN_MP, $metapartitions) || die "\nerror 663, cant open -metapartitions file:$metapartitions\n";
	print "\tfile opened.\n";
	my $in_raxml_format_section = 0;$raxml_format_partit_number = 0;
	while (my $line = <IN_MP>)
		{
		$line =~ s/\n//;$line =~ s/\r//;
# 1      | LG+G+F     | 1473       | ded04e4e9dfc8377c8ed3cd55ab61d20 | gene111_insectaNUCL_EOG55DV52_clo_pruned_fa_RMmulti, gene0_insectaNUCL_EOG50000K_clo_pruned_fa_RMmulti, gene12_insectaNUCL_EOG50K6FD_clo_pruned_fa_RMmulti
		if($line =~ /^(\d+)\s+\|\s+(\S+)\s.+\|\s+([^\|]+)$/)
			{
			my $metapart_number = $1;my $model = $2; my $partitins = $3;$partitins =~ s/\s+//g;
			$model_of_metapartition{$metapart_number} = $model;$count_metapartitions++;
		#	print "metapart_number:$metapart_number model:$model partitins:$partitins\n";
			my @split_partition = split /\,/ , $partitins;
			foreach my $part(@split_partition)
				{$partit_belongs_to_metapartit{$part} = $metapart_number};
			$meta_partitions{$metapart_number} = $partitins;
			};

		# in addiotn to partition definitions with gene IDs,
		# need to read raxml format lines, since the model codes are different.

#RaxML-style partition definitions
#Warning: RAxML allows for only a single model of rate heterogeneity in partitioned analyses. I.e. all partitions must be assigned one of three types of model: No heterogeneity (e.g. GTR); +G (e.g. GTR+G); or +I+G (e.g. GTR+I+G). If the best models for your datasetcontain different types of model for different subsets you will need to decide on the best rate heterogeneity model before you run RAxML. If you prefer to do things more rigorously, you can run separate PartitionFinder analyses for each type of rate heterogenetity Then choose the scheme with the lowest AIC/AICc/BIC score. Note that these re-runs will be quick!
#
#LGF, Subset1 = 57783-58203, 1-587, 6876-7340

		if($line =~ /^Warning. RAxML allows for only a single model/ && $in_raxml_format_section == 0)
			{$in_raxml_format_section = 1};

		if($in_raxml_format_section == 1 && $line =~ /^MrBayes/)
			{$in_raxml_format_section = 2};

		if($in_raxml_format_section == 1 && $line =~ /^\S+\,\s/)
			{
			$raxml_format_partit_number++;
			if($line =~ /^([^\,]+)\,/){$raxml_codes{$raxml_format_partit_number} = $1};
			};


		};
	close IN_MP;
	unless($count_metapartitions >= 1 ){die "\nerror, no meta-partitions read from file:$metapartitions\n"};
	print "file has been read, count metapartitions:$count_metapartitions\n";
	}else{
	};


if($arguments =~ /\-genus_level/)
	{
	$genus_level		= 1;
	}else{
	$genus_level		= 0;
	};
# $genus_level		= 0; # default 0. if 1, will trim of species namse if present, and concat at the genus level




if($arguments =~ /\-data_type\s+(\S+)/)
	{
	$force_datatype = $1;
	print "\ndata type for partition definitions forced to $force_datatype\n";
	unless($force_datatype eq "BIN" || $force_datatype eq "PROT" || $force_datatype eq "DNA")
		{die "user specified force datatype as $force_datatype, however only DNA, PROT or BIN are accepted. quitting.\n\n"};
	}else{
	print "\n-data_type not defined in your command, will try and autodetect (Warning: might be unreliable).\n";
	};


if($arguments =~ /\-required_data\s+(\S+)/)
	{
	$required_data = $1;

	my $number_of_loci = scalar @input_array;
	my $number_of_partitions = length($required_data);
	print "\nyou have specified $number_of_loci file names, and $number_of_partitions partitions in -required_data option ...\n";
	if ($number_of_loci == $number_of_partitions)
		{print "match. continuing.\n"
		}else{
		print " error, please specify same number of partitions " , 
		"(number of single alignment files should be same as the number of characters " , 
		"in the \'-required_data\' command line option).\n" , 
		"note, you can run the script without required_data, if you dont wish to filter by gene presence.\n" , 
		" quitting.\n\n";
		die "";
		};

	}else{
	$required_data = "";	
	}
# $required_data = ""; # default, nothing in this object. 
					# string like this "Y--", means you are concatenating 3 genes
					# but you only want to print a concatenate if it possesses a seq of the first gene
					# NOTE: use Y for reqruired gene for making a concatenate, 
					# otherwise use your missing data character for gene not neccessary
					# i.e. not literally '-'. 


print "

concatenate_v2.pl

settings as follows

number of matrices:" , scalar @input_array , "
missing_data_char:$missing_data_char
remove_accession:$remove_accession
genus_level:$genus_level
required_data:$required_data


";


};

##########################################################################################################################



