




# 
# 
# 
# script opens these files.
# theres no main output file as such, the result is printed to screen.
# 
# 
# 
# open(IN, $in) || die"\nerror\n";
# open(NODES , "nodes.dmp") || die "cant open nodes.dmp
# open(NAMES , "names.dmp") || die "cant open names.dmp
# open(LOG, ">$in.taxon_parse_genbank_summary_results_LOG") || die "";
# 
# 
# 
# change log 
# 2017-may-17: added option to filter at species level.
# 2018-jan-21: accession string printed to file, multiple lines where >300
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
# 
######################################################################################################




my $in 			= $ARGV[0];
my $starting_node 	= $ARGV[1];
my $level 		= $ARGV[2]; # $level = 3;# 1= order; 2=family; 3=genus; 4=species












unless($starting_node =~ /\w/){die "\nplease give starting node (NCBI taxonomy number), 
root (cellular organisms) is 131567 if you are in a rush\n"};


unless($level =~ /[1234]/)
{
die "\ncommand error\nperl script_name.pl [input] [taxonomic_starting_node] [taxonomic level]\n";
};


# from parse_ncbi_tax_db ...
$ignore_subspecies			= 1;	# 0=read taxonomies beyond species level (e.g. include subspecies, varietas, and unnamed bacteria ranks)
						# if 1, these extra ranks are ignored
						# this applies to full species name, assigned to given ncbi tax number
# globals
%ncbi_nodes;	# ncbi taxonomy structure is recorded in this object
###############
store_nodes();#
###############
##################
parse_namesfile();#
###################
$starting_name = $ncbi_nodes{$starting_node}{name};$starting_name_ff = $starting_name;
print "name for the starting node ($starting_node) which has been specified:$starting_name\n";
print "traversing NCBI taxonomy tree from this node. recording lineage for each taxon.\n\n";
$user_specified_taxon = $starting_name; $starting_name =~ s/^(\w).+$/$1/;
#################################################
traverse_nodes($starting_node , $starting_name);#
#################################################















open(IN, $in) || die"\nerror cant open file:($in)\n";

$entry=0;

while(my $line  =<IN>)
	{

#2033. Erythrobacter sp. s21-N3, complete genome
#3,012,400 bp circular DNA 
#NZ_CP011310.1 GI:890444402

#	print "$line";


	if($line =~ /^(\d+)\.\s([A-Z][a-z]+)\s([a-z]+)/)# line 1
		{
		my $entry_number = $1;#print "entry_number:$entry_number entry:$entry\n";
		my $gen = $2; my $sp = $3;
		$genera{$entry} = $gen;
		$species{$entry} = $gen . "_" . $sp;

		$store_gen{$gen . "_" . $sp}++;
		$store_sp{$gen . "_" . $sp}++;

		unless($entry_number == $entry) {print "\nerror 95 ..\n"};
		}elsif($line =~ /^[0-9\,]+\sbp/){# line 2
		
		}elsif($line =~ /^(.+)\sGI\:(.+)/)#line 3
		{
		my $accession = $1;
		$accessions{$entry}=$accession;
		}elsif($line =~ /^$/)# spacer
		{$entry++;
		}
		else{print "errror, cant parse:$line";}
	#print $line;
	}

close IN;




# lotz of these
#errror, cant parse:1415. Candidatus Tremblaya princeps PCVAL, complete genome
#errror, cant parse:1444. Mutant Legionella pneumophila subsp. pneumophila str. Hextuple_2q, complete genome



my @keys = keys %store_sp;
print "\nfound ",$#keys , " species\n";



open(LOG, ">$in.taxon_parse_genbank_summary_results_LOG") || die "";

my @ntries = keys %genera;
my @ntries2 = keys %species;

#	@keys5 = @ntries;
@keys5 = @ntries2;


foreach my $entry(@keys5)
{
my $tax = $species{$entry};
my $ac = $accessions{$entry};
$test_species = $tax;

my $test_taxonomy = $complete_lineage_for_this_species{$test_species};


unless($test_taxonomy =~ /\w+/)
	{
	# no taxonomy found for species name, try genus:
	my $genusname = $test_species; $genusname =~ s/^([A-Z][a-z]+)_[a-z]+/$1/;
	$top_taxonomy = $complete_lineage_for_this_species{$genusname};
	print "error 159 no taxonomy found for species:($test_species), try genus ($genusname)\n";
	if($top_taxonomy =~ /\w+/)
		{
		
		
		
		}else{
		print "\nerror 1258, no taxonomy found for species:($test_species), nor genus:($genusname)\n";
		#die "did you specify the correct taxonomic node? ($user_specified_taxon)\n";
		};
	};



# print "entry:$entry tax:$tax ac:$ac test_taxonomy:$test_taxonomy\n";




if($test_taxonomy =~ /\w+/)
	{

my $f1 = "NA";my $o1 = "NA";my $g1 = "NA";my $s1 = "NA";
if($test_taxonomy =~ / family:(\w+)\s/)
	{
	$f1 = $1;
	$storefams{$f1}  = $ac;
	};

	if($test_taxonomy =~ / order:(\w+)\s/)
		{
		$o1 = $1;
		$storeorders{$o1} = $ac;
		#print "entry:$entry tax:$tax lineage:$test_taxonomy\n";
		};
	if($test_taxonomy =~ / genus:(\w+)\s/)
		{
		$g1 = $1;
		$storegenus{$g1} = $ac;
		#print "entry:$entry tax:$tax lineage:$test_taxonomy\n";
		};
	if($test_taxonomy =~ / species:([A-Z][a-z]+_[a-z]+)/)
		{
		$s1 = $1;
		$storespecies{$s1} = $ac;
		#print "entry:$entry tax:$tax lineage:$test_taxonomy\n";
		};


	print LOG "entry:$entry tax:$tax ac:$ac fAM:$f1 Order:$o1 genus:$g1 species:$s1 lineage:$test_taxonomy\n";
	}





};


@fams = keys %storefams;
print scalar @fams , " fams found \n";
@genera = keys %storegenus;
print scalar @genera , " gen found \n";
@species_array = keys %storespecies; #######################################
print scalar @species_array , " species found \n";

@os = keys %storeorders;@os = sort @os;
print scalar @os , " orderss found \n";

my $accessionstring = "";


print "\n\n";
if ($level == 1)# 1= order; 2=family; 3=genus
	{
	foreach my $order(@os)
		{$ac = $storeorders{$order};
		print "order:$order ac:$ac\n";
		$accessionstring .= "$ac,"}
print "for level 1 (order), ";
	}
if ($level == 2)# 1= order; 2=family; 3=genus
	{
	foreach my $order(@fams)
		{$ac = $storefams{$order};
		print "fam:$order ac:$ac\n";
		$accessionstring .= "$ac,"}
print "for level 2 (family), ";
	}

if ($level == 3)# 1= order; 2=family; 3=genus
	{
	foreach my $order(@genera)
		{$ac = $storegenus{$order};
		print "genus:$order ac:$ac\n";
		$accessionstring .= "$ac,"}
	print "for level 3 (genus), ";
	}

if ($level == 4)# 1= order; 2=family; 3=genus
	{
	foreach my $order(@species_array)
		{
		$ac = $storespecies{$order};
	#	print "TEST1 species:$order ac:$ac\n";
		$accessionstring .= "$ac,";
		}
	print "for level 4 (species), ";
	}


print "there is one representative retained for each.\nhere is a list of their accessions:\n";


$accessionstring =~ s/\,$//;




# recently, strings are too large to be taken in one go by ncbi query box,
# awkward to manually split, so do that here and print to file.
# limit is roughly 500 , if i remember correctly

my @split_accession_string = split /\,/ , $accessionstring;

my @accession_strings_array = ();
my $current_index = 0;my $current_accession_string = "";
my $accession_strings_array_index = 0;

foreach my $accession(@split_accession_string)
	{
	$current_accession_string .= "$accession,";
	

	if($current_index >= 289)
		{
		$current_accession_string =~ s/\,$//;
		push @accession_strings_array , $current_accession_string;	
	#	 $accession_strings_array[$accession_strings_array_index] = $current_accession_string;	
	#	$accession_strings_array_index
		$current_accession_string = "";
		$current_index = 0;
		}else{
		$current_index++;
		};
	};

# last one less than 300:
$current_accession_string =~ s/\,$//;
push @accession_strings_array , $current_accession_string;	


open(OUT9 , ">tax_parse_genbank_summary_ACCESSIONS") || die "\nerror 321, \n";

foreach my $string(@accession_strings_array)
	{print OUT9 "$string\n";
	};
close OUT9;

# print "$accessionstring\n\n";

my @count_accession_srting = split /\,/ , $accessionstring;





print "
number of accession:$#count_accession_srting
node given by user:$starting_name_ff
FIN.
";
exit;

#####################################################################################################################################################
#
#
#
#####################################################################################################################################################



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
	
	
	my $name_string = $ncbi_nodes{$child}{name};$name_string =~ s/\s+/_/; ####################### sep2013
	my $child_complete_lineage = $ncbi_nodes{$current_node}{complete_lineage} . "$rank:$name_string ";#prev assigned $name
	
	$ncbi_nodes{$child}{complete_lineage} = $child_complete_lineage;
	
	$complete_lineage_for_this_species{$name_string}=$child_complete_lineage;

	$ncbi_tax_number_for_this_species{$name_string}=$child;

#	print "complete lineage:$ncbi_nodes{$child}{complete_lineage}\n";

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






