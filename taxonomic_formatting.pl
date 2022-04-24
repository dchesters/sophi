#!/usr/bin/perl



# 
# original output format for fasta file seems to be that format they use in prokaryote profilling, and SAP
# 
#
# species_to_complete_taxonomies RENAMED taxonomic_formatting.pl
#
# -seqfile XX -output XX -node XX -fasta or -newick
# 
# 
# CHANGE LOG
# 2015-09-29: does newick termnials 
# 2015-10-07: does it better
# 2015-11-03: improved format of fasta ID
# 2015-05-19: doesnt omit seqs with only partial taxonomic labelling
# 2018-01-07: prints an error message
# 2018-05-17: prints warning is user specifies trim accession for fasta file, with file that doesnt have accessions
#    		(in which case the IDs can get truncated).
# 2019-01-05: should work whether newick labeled with just genus names, or binomials.
# 2019-08-08: tiny change, another truncate string option.
# 2020-09-15: writes mothur format outfiles
# 2020-11-25: moved option $trim_accessions_from_input to command line as -trim_accessions_fasta, i frequent need to change
# 2021-01-26: if nothing found in the NCBI file for the users fasta file, this could be due to misspecified start node,
# 		thus user is reminded start taxon for node, in this case. 
# 2021-07-19: bugfix regarding $remove_branchlengths_from_newick.
# 		added useful error message for common mistake on settings.
# 2021-08-19: uses non-binomials for trees with branchlengths
# 2021-10-12: error message if user specifies newick but gives a fasta file.
# 		Now reads ID format used for MOTU in ITOL3
# 2021-11-09: wont now get stuck encountering certain weird labels
# 
# 
# 
# 
# 
# 
# 
# curiously, newick only works if adding more than one taxonomic string
# 
# 
# iv noticed my old NCBI taxonomy code (species_to_complete_taxonomies, taxonomic_report) doesnt read subgenus correctly. try instead taxon_table.pl 
# 
# 
# 
# 
# 
############################################################################################################




my $arg_string  = join ' ', @ARGV;

# additional user options:

$verbose 		= 1; # 0 simple screen output, 1=verbose, 2=ulta-verbose
$filter_level 		= 1;# 1= remove member if it is of a genus not found in tax db. 2= rm if member is of species ....

$process_backbone_tree_terminal_IDs = 0;
	# if == 1, then process ids in backbone tree, e.g. Coleoptera_Apatides_fortis to Apatides_fortis
	# if == 2, then trim off species names

# $trim_accessions_from_input = 1; # applies to fasta only
$truncate_new_taxonomic_strings = 0;# 1 = only to order, 2 === all ranks
$fasta_ID_format = 1; # applies to fasta IDs, not newick. 1= SAP, 2 = regular

					# if you have branchlengths, then choose 0 here
$remove_branchlengths_from_newick = 1; # default 1 .... 
					# ==1 even if you dont have branchlengths in input tree.
$uppercase_tax_strings = 0;

if($fasta_ID_format == 1)
	{# presumably SAP format required

	unless($filter_level == 2)
		{
		print "\nWARNING ... you have set to print entries for which lineage not found in db .... SAP will crash if these are present\n";
		};
	
	};



#@ranksarray 	= ("kingdom" , " phylum" , " order", " suborder", " superfamily",  " family", " subfamily"," tribe", " genus", " subgenus", " species group");
#my @ranksarray 	= ("kingdom" , " phylum" , " order", " suborder", " superfamily",  " family", " subfamily"," tribe", " genus", " subgenus", " species group");
#@ranksarray 	= ("kingdom" , " phylum" , " order", " family");
#@ranksarray 	= (" class" , " order");
# @ranksarray 	= (" order", " family");#, " family"

# @ranksarray 	= (" order", " suborder", " superfamily",  " family", " subfamily"," tribe", " genus", " subgenus", " species group");
# @ranksarray 	= (" order", " suborder", " infraorder", " superfamily",  " family");

# @ranksarray 	= (" class" , " order", " family", " subfamily"," tribe", " genus", " subgenus", " species group");#, " family"

# @ranksarray 	= (" order"," family");#, " family"

 @ranksarray 	= (" order", " suborder", " superfamily",  " family", " subfamily"," tribe", " genus");

# @ranksarray 	= (" subfamily"," tribe");#, " family"

# @ranksarray 	= (" order", " family", );
# @ranksarray 	= (" superfamily",  " family", " subfamily");
# @ranksarray 	= (" order", " suborder", " superfamily",  " family", " subfamily"," tribe", " genus", " subgenus", " species group");
# @ranksarray 	= (" order", " superfamily",  " family", " subfamily"," tribe");



##########################################################################################################################################


# globals:
%check_hash4;

$treefile;
$keyfile;
$fas_file;
$reference_file;
$output_filename;# 	= "$treefile.diet_clades";
%species_names;
$starting_node;
$support_cutoff = "NA";


#####################################
read_command_arguments($arg_string);#
#####################################




###############################################################################################################################################

# from parse_ncbi_tax_db ...

$ignore_subspecies			= 1;	# 0=read taxonomies beyond species level (e.g. include subspecies, varietas, and unnamed bacteria ranks)
						# if 1, these extra ranks are ignored
						# this applies to full species name, assigned to given ncbi tax number
# globals

%ncbi_nodes;	# ncbi taxonomy structure is recorded in this object

%assign_these_new_terminals_to_taxon;


###############
store_nodes();#
###############



###################
parse_namesfile();#
###################


$starting_name = $ncbi_nodes{$starting_node}{name};$starting_name_ff = $starting_name;
print "name for the starting node ($starting_node) which has been specified:$starting_name\n";
print "traversing NCBI taxonomy tree from this node. recording lineage for each taxon.\n\n";
$starting_name =~ s/^(\w).+$/$1/;

#################################################
traverse_nodes($starting_node , $starting_name);#
#################################################



#		my $thing = $complete_lineage_for_this_species{$name};
#		#print "name:$name thing:$thing\n";
#		if($thing =~ /order\:$order\s.+\sspecies\:([A-Z][a-z]+_[a-z]+)\s/)	



###############################################################################################################################################



my %terminals = ();
my $root_identity = "";
my $date = localtime time;

%query_IDs;	# the fasta file contains names of new entries to be placed onto backbone tree. 
		# each are stored in the hash %query_IDs.


if($assign_to_tree==1)
	{
	###############	
	read_newick();#
	###############	

	if($this_file_is_fasta_not_newick >= 10)
		{
		print "\n\n\nERROR, you specified newick input file but pretty sure its fasta. quitting.\n";
		die "\n\n";
		};

	my @ranks3 = keys %rank_assignment_counts;
	foreach my $rank(@ranksarray)
		{
		my $no = $rank_assignment_counts{$rank};print "rank:$rank number assigned:$no\n";	
		};

	print "user speified base node for taxonomies:$starting_name\n";
	print "\n\n\nnewick read and taxonomies assigned. end of script\n";
	exit;
	
	}else{

	############
	read_fas();#
	############

	print "user speified base node for taxonomies:$starting_name\n";
	print "\n\n\nFASTA read and taxonomies assigned. end of script\n";
	exit;

	};







#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################






sub read_fas
{

open(IN, $fas_file) || die "\nerror ... fasta file you have given ($fas_file) cannot be opened. 
make sure it is named correctly and in the working directory. quitting\n";

my $count_diet_IDs = 0;
open(OUT1, ">$output_filename") || die"\nerror 140\n";
open(MOTHURREF, ">mothur.reference") || die"\nerror 221\n";
open(MOTHURTAX, ">mothur.taxonomy") || die"\nerror 222\n";
open(METACODER, ">metacoder_format") || die"\nerror 246\n";
open(METACODER_GENUS, ">metacoder_genus") || die"\nerror 252\n";


my $fasta_entry = "";
while(my $fasta_line= <IN>)
	{
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


unless(length($fasta_entry)<=2)
	{
	$fasta_entry =~ s/^>//;
	########################################
	process_entry($fasta_entry);#
	########################################
	}

close IN;
close OUT1;
close MOTHURREF;
close MOTHURTAX;
close METACODER;
close METACODER_GENUS;

sub process_entry
{
my $line = shift;
my $id = "";
if ($line =~ /^(.+)\n/ )
	{
	$id = $1;
	$line =~ s/^.+\n//;#$line =~ s/\n//;$line =~ s/\r//;my $current_seq_length = length($line);
	}else{
	die "\nerror 5631.\n"
	};

if($trim_accessions_from_input == 1)
	{
	if($id =~ /^[A-Z][a-z]+_[a-z]+$/)
		{
		$warning_count_1++;
		if($warning_count_1 < 100)
			{
			print "warning, you specified trim accessions, but this member seems to have none, might truncate:$id\n"
			}elsif($warning_count_1 == 100)
			{
			print " ...... warning print limit reached ...... \n"
			};
		};

	$id =~ s/^([A-Z][a-z]+_[a-z]+).+/$1/;
	};

my $genus = $id;#print "id:$id\n";
		if($genus =~ s/^([A-Z][a-z]+)_[a-z]+$/$1/)
			{
			$count_binomials++;
			}elsif($genus =~ s/^[a-z]+_([A-Z][a-z]+)_[a-z]+/$1/)
				{
				$extended_ids++;
				}elsif($genus =~ s/^([A-Z][a-z]+)_[a-z]+_[0-9]+$/$1/){
				#Lasiopogon_cinctus_158267666
				}
				elsif($genus =~ s/^([A-Z][a-z]+)_[a-z]+_.+/$1/)
				{
				# Xylocopa_virginica_EU271670
				}
				elsif($genus =~ s/^([A-Z][a-z]+)_SQ[0-9]+$/$1/)
				{
				# format used for motu in ITOL3: Zythos_SQ39849
				}else{
				unless($filter_level == 0){print  "\nWARNING 153. not parsed iD:$genus\n"};
				};

		$count_diet_IDs++;
		$query_IDs{$id}= 1;

		my $look_for_ID; # 1= remove member if it is of a genus not found in tax db. 2= rm if member is of species ....
		if($filter_level == 0){$look_for_ID = $genus};
		if($filter_level == 1){$look_for_ID = $genus};
		if($filter_level == 2){$look_for_ID = $id};

		if(exists($complete_lineage_for_this_species{$look_for_ID}))
			{
			my $completelineage = $complete_lineage_for_this_species{$look_for_ID};
		#	print "completelineage:$completelineage\n";

			my $phylum; my $order; my $fam; my $subfam; my $tribe; my $genus; my $species;

			my @lineage_array;my @lineage_array2;# superkingdom was not parsed if starting \s
			my @lineage_array3;
			foreach my $current_rank (@ranksarray)
				{
				if($completelineage =~ /($current_rank):(\S+)/i)
					{
					my $rank = $1; my $tax = $2; $tax =~ s/_/ /g;
					my $newstringX = $rank . ": " . $tax;push @lineage_array , $newstringX;

				#	if($truncate_new_taxonomic_strings == 1)
				#		{$tax =~ s/^(\w\w\w\w\w\w)\w+/$1/};
					if($truncate_new_taxonomic_strings == 1)
						{
					#	$tax =~ s/^(\w\w\w\w\w\w\w\w)\w+/$1/;
						if($rank eq "order")
							{
							print "truncating $tax\n";
							$tax =~ s/^(\w\w\w)\w+/$1/;
							$tax = uc($tax);
							};
						}elsif($truncate_new_taxonomic_strings == 2)
						{
							print "truncating $tax\n";
							$tax =~ s/^(\w\w\w\w\w\w)\w+/$1/;
							# $tax = uc($tax);
						};

					my $newstringY = $tax . "_";push @lineage_array2 , $newstringY;
					my $newstringZ = $rank . "__" . $tax . ";";push @lineage_array3 , $newstringZ;
# 114	k__Archaea;p__Euryarchaeota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter;


					};
				}
			my $linstring = join ', ', @lineage_array;my $linstring2 = join '', @lineage_array2;my $linstring3 = join '', @lineage_array3;
			my $linstring4 = $linstring3;$linstring4 =~ s/\s+//g;
#			print "id:$id, lineage:$linstring\n";

			# presumably underscores were removed here to conform to sap format, so add back in afterwards
			my $binomial = $id; $binomial =~ s/_/ /g;#print "binomial:$binomial\n";

			$counter65++;my $abritratry = "DC$counter65";

			my $sap_format 		= "$abritratry ; " . $linstring . ", species: $binomial ; $binomial";
			my $binomial3 		= $binomial;$binomial3 =~ s/\s+/_/;	
			my $mothur_format 	= "$abritratry\t$linstring3" . "species__$binomial3;";
			my $metacoder_format 	= "$abritratry\tclass__Insecta;$linstring4" . "species__$binomial3;";
			my $metacoder_format2 	= "$abritratry\tclass__Insecta;$linstring4";


			#print "\nsap_format:$sap_format\n";
			$binomial =~ s/ /_/g;

#>mySequenceID ; family: Hominidae, genus: Homo, species: Homo sapiens ; Myself
#The following taxonomic levels are recognized: kingdom, phylum, class, order, suborder, superfamily,  family, subfamily, supertribe, tribe, subtribe, supergenus, genus, subgenus, species group, species subgroup, species, subspecies. The name after the last semicolon is the organism name. This is sometimes different from the species or subsp

			if($fasta_ID_format == 1)
				{
				# print "sap_format:$sap_format\n"; # Cubitalia_morio
				if($sap_format =~ /\s[a-z]+\:/) # $ranksarray[0]
					{
					print OUT1 ">$sap_format\n$line";

					print MOTHURREF ">$abritratry\n$line";
					print MOTHURTAX "$mothur_format\n";
					print METACODER "$metacoder_format\n";

					if($printedalready{$linstring4} == 1)
						{
						
						}else{
						print METACODER_GENUS "$metacoder_format2\n";
						};
					$printedalready{$linstring4} = 1;

				# 13	k__Archaea;p__Euryarchaeota;c__Thermoplasmata;o__E2;

					}else{
					print "Warning 355: $sap_format\n"
					};
				};
			if($fasta_ID_format == 2)
				{
				print OUT1 ">$linstring2" , "_" , $binomial , "\n$line";
					
				};
			$present++;#print "fasta ID is PRESENT from taxonomy DB:($id)\n";
			}else{
			if($filter_level == 0)
				{
				print OUT1 ">$id\n$line";
				}else{
				$absent++; print "filter_level = $filter_level, fasta ID is MISSING from taxonomy DB:($id), look_for_ID:($look_for_ID)\n";
				};
			}

};



print "
fasta IDs present in tax DB:$present
fasta IDs absent  in tax DB:$absent
";

print "
count_binomials:$count_binomials 
extended_ids:$extended_ids
";

unless($present >= 2)
	{
	print "\nerror, nothing found. did you specify starting node correctly? node:$starting_node is taxon:$starting_name_ff\n\n";
	};


}#sub read_fas





#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################








#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


sub get_species_list
{
my $list_of_terminals = shift;

my %hash_of_terminals = ();

if($verbose == 1){print "list_of_terminals:$list_of_terminals\n"};

my @array_of_terminals = split /\t/, $list_of_terminals;

foreach my $termnal(@array_of_terminals)
	{
	$termnal =~ s/^([^_]+)_.+/$1/;
	$hash_of_terminals{$termnal}= 1;
	}

my @array_of_terminals2 = keys %hash_of_terminals;

if($verbose == 1){
#print "list_of_species:@array_of_terminals2\n";
}

my $list_of_species = join ' ', @array_of_terminals2;

return($list_of_species);

}



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

# -seqfile -output -node -fasta

if($arguments =~ /-seqfile\s+(\S+)/)
	{
	$fas_file = $1;
	}else{
	}

if($arguments =~ /-trim_accessions_fasta\s+(\S+)/)
	{
	$trim_accessions_from_input = $1; # applies to fasta only
	}else{

	if($arguments =~ /-fasta/)
		{
		die "\n	command error, missing -trim_accessions_fasta [moved to command line 2020-11]\n\n	";
		}else{
		$trim_accessions_from_input = 0;	
		};
	}



if($arguments =~ /-fasta/)
	{
	$inputformat = "fas";
	}else{
	if($arguments =~ /-newick/)
		{
		$assign_to_tree=1;
		}else{
		die "error";
		};
	}


if($arguments =~ /-output\s+(\S+)/)
	{
	$output_filename = $1;
	}else{
	}



if($arguments =~ /-node\s+(\d+)/)
	{
	$starting_node = $1;
	}else{
	print "user did not given NCBI taxonomic number. using default 33208 (Metazoa)
if you have references outside this default, you need to get the appropriate NCBI taxonomy number from http://www.ncbi.nlm.nih.gov/taxonomy/
and input this here with the option -node taxon_number
";
	$starting_node = 33208;
	}


#$output_filename	= "$treefile.query_clades";


print "\n
user options have been read.....
support_cutoff:			$support_cutoff
taxon node encompassing refs:	$starting_node
treefile:			$treefile
query fasta file:		$fas_file
reference_file:			$reference_file

output will be written to file:	$output_filename


";



}#sub read_command_arguments



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



sub read_keyfile
{
open(IN, $keyfile) || die "\n\nerror 91. cant open $keyfile\n";
my $species_strings_read=0;

while (my $line = <IN>)
	{
	#print "\n$line";

	#MeeRBBegrac Berberis_gracilis 258166 species no_rank:eudicotyledons no_rank:eudicotyledons order:Ranunculales family:Berberidaceae genus:Berberis species:gracilis

	if($line =~ /^(\S+)\s(\S+)\s(\S+)\s(\S+)\s(.+)/)
		{
		my $tobycode = $1;my $species_name = $2;my $ncbi_number = $3; my $rank = $4; my $taxonomic_path = $5;
		#print "\ttobycode:$tobycode species_name:$species_name ncbi_number:$ncbi_number rank:$rank taxonomic_path:$taxonomic_path\n";
		$species_names{$tobycode} = $species_name;$taxonomic_paths{$tobycode} = $taxonomic_path;
		$species_strings_read++;
		$tobycodes{$species_name} = $tobycode;

	#	if($tobycode eq "MLAA4"){print "$tobycode $species_name ... quit\n\n";die ""}

		}
	}

close(IN);

if($species_strings_read == 0){die "\n\nerror 112. seems to be problem reading keyfile, lines in that file dont match regex.\n"}

print "$species_strings_read species strings read from file, plus associated taxonomic information\n";

}




#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



sub store_nodes
{

my %rank_hash;

# nodes.dmp contains each taxon, described in a tree-like structure. 
# so each node has a name (eg polyphaga) and the higher group node to which it belongs (beetles). 

open(NODES , "nodes.dmp") || die "cant open nodes.dmp
go to ftp://ftp.ncbi.nih.gov/pub/taxonomy/\nget the file taxdump.tar.gz\nunzip and place in working directory\nquitting.\n";

print "\nreading files from NCBI taxonomy database .... nodes.dmp .... ";


my $line_counter=0;
while (my $line = <NODES>)
	{
	if($line =~ /^(\d+)\t[\|]\t([^\t]+)\t[\|]\t([^\t]+)\t[\|]\t/)
		{
		my $tax_id = $1;my $parent_tax_id = $2;my $rank = $3;
		$rank_hash{$rank}++;
#		print "tax_id:$tax_id parent_tax_id:$parent_tax_id rank:$rank\n";

			$ncbi_nodes{$tax_id}{rank} = $rank;
		#	$ncbi_nodes{$tax_id}{rank_code} = $current_rankcode;
			$ncbi_nodes{$tax_id}{parent} = $parent_tax_id;
			$ncbi_nodes{$parent_tax_id}{child_nodes} .= "\t" . $tax_id;

		}else{
		print "line_counter:$line_counter line:$line";
		die "UNEXPECTED LINE:$line\nquitting\n";
		}
	$line_counter++;
	}

close(NODES);

#my @ranks = keys %rank_hash;@ranks = sort @ranks;
#print "ranks found in nodes.dmp:\n";
#print LOG "ranks found in nodes.dmp:\n";

#foreach my $rank(@ranks){print "$rank\t" , $rank_hash{$rank} , "\n";print LOG "$rank\t" , $rank_hash{$rank} , "\n"};

my @all_nodes = keys %ncbi_nodes;@all_nodes = sort @all_nodes;

print scalar @all_nodes , " nodes have been read.\n";


}




#####################################################################################################
#
#
#
#####################################################################################################


sub parse_namesfile
{

# here just parse the scientific name of each node. ignore synonyms etc

open(NAMES , "names.dmp") || die "cant open names.dmp
go to ftp://ftp.ncbi.nih.gov/pub/taxonomy/\nget the file taxdump.tar.gz\nunzip and place in working directory\nquitting.\n";


print "\nnames.dmp, parsing 'scientific name', ignoring others ... ";

my $names_line_counter=0;
while (my $line = <NAMES>)
	{
# 24	|	Shewanella putrefaciens	|		|	scientific name	|

	if($line =~ /^(\d+)\t[\|]\t([^\t]+)\t[\|]\t([^\t]*)\t[\|]\tscientific name/)
		{
		my $tax_id = $1;my $name = $2;#my $rank = $3;
		# print "tax_id:$tax_id name:$name\n";

		# if you want to remove non-alphanumerical characters from assigned species names:
		$name =~ s/[\(\)\,\[\]\'\#\&\/\:\.\-]/ /g;
		$name =~ s/\s\s+/ /g;$name =~ s/\s+$//;$name =~ s/^\s+//;
		$ncbi_nodes{$tax_id}{name} = $name;

		$names_line_counter++;#print "$names_line_counter\n";


		}else{
		if($line =~ /^(\d+).+scientific name/){die "UNEXPECTED LINE:\n$line\nquitting\n"}
		}

	}

close(NAMES);

print "$names_line_counter names parsed.\n";




}



#####################################################################################################
#
#
#
#####################################################################################################




sub traverse_nodes
{
my $current_node = $_[0];my $current_node_taxstring = $_[1];# $current_node_taxstring to be deprecated

my $child_nodes = $ncbi_nodes{$current_node}{child_nodes};$child_nodes =~ s/^\t//;
my @child_nodes_array = split(/\t/, $child_nodes);



if($current_node == $starting_node)
	{
	my $rank = $ncbi_nodes{$starting_node}{rank};$rank =~ s/\s/_/g;$how_many{$rank}++;
	my $current_node_taxstring = substr $current_node_taxstring, 0,1;
	my $originalname = $ncbi_nodes{$starting_node}{name};$originalname =~ s/[\s\t]/_/g;
	my $child_complete_lineage = $ncbi_nodes{$starting_node}{complete_lineage} . "$rank:$originalname ";
	}



foreach my $child(@child_nodes_array)
	{
	unless($current_node == $child)# one the child nodes of the root node (1), is also 1 
	{
	my $rank = $ncbi_nodes{$child}{rank};$rank =~ s/\s/_/g;$how_many{$rank}++;
	
	
	my $name_string = $ncbi_nodes{$child}{name};$name_string =~ s/\s+/_/g; ####################### sep2013 .. fixed mar2015
	my $child_complete_lineage = $ncbi_nodes{$current_node}{complete_lineage} . "$rank:$name_string ";#prev assigned $name
	
	$ncbi_nodes{$child}{complete_lineage} = $child_complete_lineage;

	#print "storing complete lineage for taxon:$name_string\n";	
	$complete_lineage_for_this_species{$name_string}=$child_complete_lineage;

	$ncbi_tax_number_for_this_species{$name_string}=$child;

	if($name_string =~ /Zorochros/)
		{
		#print "name_string:($name_string) child:$child complete lineage:$ncbi_nodes{$child}{complete_lineage}\n";
		};


	my $originalname = $ncbi_nodes{$child}{name};$originalname =~ s/[\s\t]/_/g;

	my $name_assignment_to_taxnumber = "";

	if($ignore_subspecies == 1)
		{
		if ( $rank eq "subspecies" || $rank eq "varietas"   || $nodes{$current_node}{complete_lineage} =~ / species:/)# 
			{
			my $parentoriginalname = $ncbi_nodes{$current_node}{name};
			if($parentoriginalname =~ /^[A-Z][a-z]+\s[a-z]+$/)
				{
				$parentoriginalname =~ s/[\s\t]/_/g;
				$name_assignment_to_taxnumber = $parentoriginalname;
			#	print "node:$child named:$originalname rank:$rank appears to be subspecfic\n";
			#	print "\tassigning parent name instead:$parentoriginalname\n\n";
				}else{$name_assignment_to_taxnumber = $originalname}
			}else{
			$name_assignment_to_taxnumber = $originalname
			}

		}else{
		$name_assignment_to_taxnumber = $originalname
		}

	#print "$name_assignment_to_taxnumber $child $rank $child_complete_lineage\n";


		###########################################
		traverse_nodes($child , $child_taxstring);#
		###########################################
	}}


	
}#sub traverse_nodes





#####################################################################################################
#
#
#
#####################################################################################################

















sub record_tree3
{
my $t5 = shift;

my $tree1= "";
open(IN, $t5) || die "\n\nerror 1408 cant open file $treefile\n\n";
while (my $line = <IN>)
	{
	if($line =~ /(.+\(.+)/){$tree1 = $1}
	}
close(IN);

print "\nlooking at tree:$treefile\n";

if($tree1 =~ s/\:[0-9\.]+//g){print "\nwarning, branhc lengths removed\n"};
$tree1 =~ s/ //g;


my $newick_string 	= $tree1;
my $interal_node	= 0;



# new newick parser ....

while ($newick_string =~ s/\(([^\(\)]+)\)/INTERNAL_NODE_$interal_node/) # rooted bipart raxml tree, the 2 last nodes will have no boot values
	{
	my $node = $1;my $boot = $2; my $nodeID = "INTERNAL_NODE_$interal_node"; #print "nodeID:$nodeID node:$node\n";

	# found a pair of parentheses (with no parentheses inside).
	# this string is taken, and replaced with a single node (INTERNAL_NODE_\d)
	# node is numbered (\d in variable $interal_node), 
	# the number being incremented at the end of the loop
	# find what is at the current node (decendents), by splitting at the commas


	my @child_nodes = split /\,/ , $node;
	$child_counts3{$nodeID} = $#child_nodes;



	for $i(0 .. $#child_nodes)
		{
		#print "$child_nodes[$i]\n";

		# record each decendent of the current node, in a hash nodes{ the current node }{ child number }

		$nodes3{$nodeID}{$i} 			= $child_nodes[$i];
		$nodes3{$child_nodes[$i]}{parent} 	= $nodeID;

		unless($child_nodes[$i] =~ /INTERNAL_NODE_/){$terminals{$child_nodes[$i]} =1}

		}
	#print "node:$interal_node\n\tchild nodes:@child_nodes\n";

	$root_node3 = $nodeID; # this is a global and will keep the value in the final loop, so the root can be identified later.

	$interal_node++;

	print "\nnewick_string:$newick_string\n";
$newick_string =~ s/\((INTERNAL_NODE_\d+)\)/$1/;

	}


print "your tree has been read, it has $interal_node nodes.\n";
#print "newick_string:$newick_string\n";die;


unless($interal_node >= 2){die "\nerror reading your phylogeny.\n"}


}#sub record_tree2




#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


sub read_newick
{

open(IN, $fas_file) || die "\nerror ... fasta file you have given ($fas_file) cannot be opened. 
make sure it is named correctly and in the working directory. quitting\n";

my $count_diet_IDs = 0;
open(OUT1, ">$output_filename") || die"\nerror 140\n";
my $fasta_entry = "";
while(my $newick_line= <IN>)
	{
	if($newick_line =~ /\(/)
		{
		my $tr = proccess_tree($newick_line);
		print OUT1 $tr;
		}else{
		print OUT1 $newick_line;
		}

	if($newick_line =~ /^>/){$this_file_is_fasta_not_newick++};

	};


####################################################################################################

sub proccess_tree
{
my $tree = shift;



# print "tree:$tree\n";

my $lineages_assigned =0;


if($remove_branchlengths_from_newick == 1)
	{
	if($tree =~ s/\:[0-9\.]+//g){print "\nwarning, branhc lengths removed\n"};

##############################################################################################################


# tree with binomials:
# while($tree =~ /[\(\)\,]([A-Z][a-z]+_[a-z]+)[\(\)\,]/)

# hack to work on tree with just genus names:
# while($tree =~ /[\(\)\,]([A-Z][a-z]+)[\(\)\,]/)

# hope this works on both:
 while($tree =~ /[\(\)\,]([A-Z][a-z]+[_a-z]*)[\(\)\,]/)

# accessions
# while($tree =~ /[\(\)\,]([A-Z][a-z]+[_a-zA-Z0-9]*)[\(\)\,]/)

	{
	my $binomial  =$1;$binomials_found_in_newick++; #print "\n";
	if($binomials_found_in_newick =~ /000$/){
		print "line 1000, $binomial. proccesed $binomials_found_in_newick binomials in newick. Binoms_which_tax_not_found:$newick_binomials_for_which_taxonomies_not_found\n";
		};
	#  if you get stuck in this loop, try adding ranks to: @ranksarray 	= (" order", " family", " subfamily");#, " family"


	my $add_tax_string = "NA_";
	my $look_for_ID = $binomial;

	my $completelineage;my $genus;
	if(exists($complete_lineage_for_this_species{$look_for_ID}))
		{
		$completelineage = $complete_lineage_for_this_species{$look_for_ID};
		}else{
		$genus = $look_for_ID;$genus =~ s/^([A-Z][a-z]+)_.+/$1/; # print "no sucess looking for binomial:$binomial, trying genus name:$genus\n";
		$completelineage = $complete_lineage_for_this_species{$genus};
		}

	if($completelineage =~ /\w/)
		{
		$add_tax_string = "";
		foreach my $current_rank (@ranksarray)
			{
			if($completelineage =~ /($current_rank):(\S+)/i)
				{
				my $rank = $1; my $tax = $2; $tax =~ s/_/ /g;
				$rank_assignment_counts{$rank}++;
			#	print "rank:$rank\n";
				if($truncate_new_taxonomic_strings == 1 && $rank =~ /order/)
					{
					$tax =~ s/^(\w\w\w)\w+/$1/;
					};

				if($truncate_new_taxonomic_strings == 2)
					{
					$tax =~ s/^(\w\w\w\w\w\w\w\w)\w+/$1/;
					};

				if($uppercase_tax_strings == 1 && $rank =~ /order/)
					{$tax  =uc($tax)};

				$add_tax_string .= "$tax" . "_";
				};
			}
		#$add_tax_string =~ s/_$//;
		$lineages_assigned++;
		unless($add_tax_string =~ /\w/){$add_tax_string = "NA_"}; # script got stuck for label Collembola_here
		}else{
		$newick_binomials_for_which_taxonomies_not_found++;	
		print "no sucess looking for binomial:$binomial nor genus name:$genus\n";

		}

	if($process_backbone_tree_terminal_IDs == 2) # if == 2, then trim off species names
		{
		$tree =~ s/([\(\)\,])([A-Z][a-z]+)_[a-z]+([\(\)\,])/$1$add_tax_string$2$3/;
		}else{

		# STRICT BINOMIALS:
	#	$tree =~ s/([\(\)\,])([A-Z][a-z]+_[a-z]+)([\(\)\,])/$1$add_tax_string$2$3/;

		# either binomial or just genus
		$tree =~ s/([\(\)\,])([A-Z][a-z]+[_a-z]*)([\(\)\,])/$1$add_tax_string$2$3/;


# hack part 2:
#		$tree =~ s/([\(\)\,])([A-Z][a-z]+)([\(\)\,])/$1$add_tax_string$2$3/;
		

# 	$tree =~ s/[\(\)\,]([A-Z][a-z]+[_a-zA-Z0-9]*)[\(\)\,]/$1$add_tax_string$2$3/;


		};

	print "";
	};

print "
binomials_found_in_newick:$binomials_found_in_newick
of these, lineages_assigned:$lineages_assigned
newick_binomials_for_which_taxonomies_not_found:$newick_binomials_for_which_taxonomies_not_found
";

if($binomials_found_in_newick >= 1)
	{
	print "binomias found\n";
	}else{
	print "Warning 1087. no binomials found \n";
	if($remove_branchlengths_from_newick == 0 && $tree =~ /\:/)
		{
	#	print "\tmaybe cause your tree has branchlengths, and you set script not to remove these?\n";
		};
	};


unless($lineages_assigned >= 1)
	{
	print "seems error, no lineages found for binomials of your tree\n";
	};





##############################################################################################################


}else{ # if($remove_branchlengths_from_newick == 1)

# here if $remove_branchlengths_from_newick == 0

if($tree =~ /\:/)
	{
	print "\$remove_branchlengths_from_newick:$remove_branchlengths_from_newick\n";
#	die "\nERROR. your tree HAS branchlengths, quitting.\n"

	}else{
	};


while( $tree =~ /[\(\)\,]([A-Z][a-z]+_[a-z]+)[\:]/ ) # this doesnt work with partial labelleds
#while( $tree =~ /[\(\)\,]([A-Z][a-z]+_[^\:]+)[\:]/ ) # so try this
	{
	my $binomial  =$1;$binomials_found_in_newick++; # print "$binomial\n";
	if($binomials_found_in_newick =~ /000$/)
		{
		print "line 1112. proccesed $binomials_found_in_newick binomials in newick. sucesful:$lineages_assigned tax_not_found:$newick_binomials_for_which_taxonomies_not_found\n";

		};

	my $add_tax_string = "NA_";
	my $look_for_ID = $binomial;

	my $completelineage;my $genus;
	if(exists($complete_lineage_for_this_species{$look_for_ID}))
		{
		$completelineage = $complete_lineage_for_this_species{$look_for_ID};
		}else{
		$genus = $look_for_ID;$genus =~ s/^([A-Z][a-z]+)_.+/$1/;
		$completelineage = $complete_lineage_for_this_species{$genus};
		}

	if($completelineage =~ /\w/)
		{
		$add_tax_string = "";
		foreach my $current_rank (@ranksarray)
			{
			if($completelineage =~ /($current_rank):(\S+)/i)
				{
				my $rank = $1; my $tax = $2; $tax =~ s/_/ /g;
				$rank_assignment_counts{$rank}++;

			#	if($truncate_new_taxonomic_strings == 1)
			#		{$tax =~ s/^(\w\w\w\w\w\w)\w+/$1/};
				if($truncate_new_taxonomic_strings == 1)
					{$tax =~ s/^(\w\w\w\w\w\w\w\w)\w+/$1/};

				$tax  =uc($tax);
				$add_tax_string .= "$tax" . "_";
				};
			}
		#$add_tax_string =~ s/_$//;
		$lineages_assigned++;
		}else{
		$newick_binomials_for_which_taxonomies_not_found++;	
		print "no sucess looking for binomial:$binomial nor genus name:$genus\n";

		}

	if($process_backbone_tree_terminal_IDs == 2) # if == 2, then trim off species names
		{
		$tree =~ s/([\(\)\,])([A-Z][a-z]+)_[a-z]+([\:])/$1$add_tax_string$2$3/;
		}else{
		$tree =~ s/([\(\)\,])([A-Z][a-z]+_[a-z]+)([\:])/$1$add_tax_string$2$3/;
	#	$tree =~ s/([\(\)\,])([A-Z][a-z]+_[^\:]+)([\:])/$1$add_tax_string$2$3/;

	
		}


	}

print "
binomials_found_in_newick:$binomials_found_in_newick
of these, lineages_assigned:$lineages_assigned
newick_binomials_for_which_taxonomies_not_found:$newick_binomials_for_which_taxonomies_not_found
";

if($binomials_found_in_newick >= 1)
	{
	print "binomias found\n";
	}else{

	print "$tree\n";

	print "Warning 1193. no binomials found \n";
	print " tip, this might happen if you have no branchlengths in input tree and you didnt set script to remove branchlengths\n";
	print " (ie script is then assuming branchlengths when looking for terminal IDs)\n\n";
	if($remove_branchlengths_from_newick == 1) 
		{
			unless( $tree =~ /\:/)
			{
		#	die "\tmaybe cause your tree has NO branchlengths, and you set script to expect these?\n";
			};
		}
	};



unless($lineages_assigned >= 1)
	{print "seems error, no lineages found for binomials of your tree\n"}




};

##############################################################################################################







return($tree);

};#sub proccess_tree

####################################################################################################














};#sub read_newick

#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




