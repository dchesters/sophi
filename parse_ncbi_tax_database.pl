

###################################################################################################################################
#
#
#
# 	parse_ncbi_tax_database.pl, 
#	Perl script for parsing taxonomic information from the NCBI taxonomy database files names.dmp + nodes.dmp
#
#    	Copyright (C) 2013-2018 Douglas Chesters
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
#	NOTES:
#	broadly similar to build_taxkeys.pl (James Abbott and Toby Hunt), 
#	but simpler, and requires no perl modules or database software.
#	
#	first reads nodes.dmp, then names.dmp
#	then traverses nodes from specified start, 
#	strings are printed to file during the process, to file such as: key_Jan2010_Coleoptera
#	taxonomic strings are described in Hunt & Vogler 2008 MPE 47:289-301
#	taxon counts printed to log
#	
#	scientific name is used, misspelling , synonym, equivalent name, common name, ignored
#
#
#	CHANGE LOG:
#  	nov2009: given names are by the filtered method (used to generate strings) rather than raw
#  	jan2010: option to ignore subsp
#  	aug2011: rm colon from BOLD id's
#  	apr2013: looks for plant subspecific rank of varietas (are there any other subspecific ranks?)
#  	aug2013: option to make keyfiles for paraphyletic groups.
#		assignment of species name to subspecies groups, apllied to binomials in addition to tobycodes
#  	sep2013: give complete name for internal nodes >1 word
#		option to build strings using internal numbers for same-named sister taxa 
#		instead of increasing substring for these
#		should be better for some downstream string comparisons, since it better differentiates taxon strings
#	oct2013: cleaned up a bit, include GPL
#	apr2016: options given in command instead of written into script, more appropriate as part of a pipeline.
#	08aug2017: minor change, couple more non-standard characters discovered in species IDs, now removed	
#	17jan2018: some more detailed error messages
#
#
#
#
#
#
#
#
###################################################################################################################################



my $arg_string  = join ' ', @ARGV;

#####################################
read_command_arguments($arg_string);#
#####################################


# $starting_node = $ARGV[0];

			# Sample taxonomy IDS:
			# 3398 = Magnoliophyta
			# Eukaryota = 	2759	, Insecta = 50557, Endopterygota 	= 33392, 	Coleoptera = 7041,
			# myxophaga = 	63907	, Lucanus = 41108, Neuroptera 		= 7516, 	caraboidea = 535382
			# diptera = 	7147	, weevils = 71529, Leps 		= 7088, 	Carabidae  = 41073
			# Cucujoidea = 71526, Chalcidoidea = 7422, Scarabaeoidea = 75546 , 
			# Papilionoidea = 37572, Apoidea = 34735, apidae=  7458. Ichneumonoidea = 7401
			# 1 = root. 2= bacteria, 2157=archaea, 131567 = cellular organisms, 4751=fungi, 33090=viridiplantae, 33208=metazoa, 
			# eukaryota  = 2759,  hymenoptera = 7399


@omit_taxa	= ();			# by default this array should contain nothing. 
					# you can make a paraphyletic group if you wish,
					# by putting ncbi taxon id's for those to be ommited, in this array
					# e.g. 'protists': euk (2759), ommiting (4751, 33090, 33208)


$date = localtime time;
$month = $date;
$month=~ s/\w+\s+(\w+)\s+\d+\s+\S+\s+(\d+)$/$1$2/;
$month= "Oct2013";
print "month:$month\n";



$ignore_atypical_ranks_to_shorten_codes = 1;	# 0==no. 1==yes, which uses the standard ranks genus, tribe, subfamily, family, superfamily, order 


# $ignore_subspecies			= 1;	# 0=build codes beyond species level (e.g. include subspecies, varietas, and unnamed bacteria ranks)
						# if 1, these extra ranks are ignored, and species strings are assigned to them 
						# this applies to both tobycode and full species name, assigned to given ncbi tax number


# $parse_species_only			= 1;	# default = 0


$rank_numbers_in_taxcodes		= 0;	# this puts numbers into the codes just to give the rank, 
						# for example ...Co7.., you know ...Co must be taxonomic string up to genus level
						# this affects downstream analyses, e.g. 2 identical seqs are more likely to also have identical ID
						# 

						 
$duplicate_taxonstring_numbers		= 1;	# 
						# two sister taxon names starting with same letter, 
						# here can use numbers to differentiate rather than increasing the substring
						# doesnt apply to species level

# globals
@all_nodes;
%nodes;
%all_taxstrings;
%code_counts;
%nonbinomial_recorder;
%binomial_recorder;
$currently_species_level_or_lower =0;
%unidentified_sp_hash;



open(LOG , ">parse_ncbi_tax_database_LOG") || die "cant open LOG\n";
print LOG "\nrunning script parse_ncbi_tax_database.pl\n$date\nstarting_node:$starting_node ignore_atypical_ranks_to_shorten_codes:$ignore_atypical_ranks_to_shorten_codes\n";
print "\n\nrunning script parse_ncbi_tax_database.pl\n$date\nstarting_node:$starting_node\n";


print "\nstoring nodes\n";

###############
store_nodes();#
###############


print "\nparse_namesfile\n";

###################
parse_namesfile();#
###################



$starting_name = $nodes{$starting_node}{name};$starting_name_ff = $starting_name;

$starting_name =~ s/^(\w).+$/$1/;

$key_file = "key_" . $month . "_$starting_name_ff";
open(OUT ,">$key_file") || die "cant open \n";
$key_file2 = "key_" . $month;
# open(OUT2 ,">$key_file2") || die "cant open \n";

print "\ntraverse_nodes\n";

#################################################
traverse_nodes($starting_node , $starting_name);#
#################################################

close(OUT);
close(OUT2);


@keys = keys %how_many;@keys = sort @keys;

print LOG "\nfrom node:$starting_name_ff\n";
print "\nfrom node:$starting_name_ff\n";

foreach $key(@keys)
	{
	print LOG "$key\t" , $how_many{$key} , "\n";
	# print "$key\t" , $how_many{$key} , "\n";
	}



print "\nend of script\n";
print LOG "\nend of script\n";

close(LOG);

my @keys = keys %binomial_recorder;@keys = sort @keys;

$count_binomials = 0;
$count_binomials2 = 0;

foreach $key(@keys)
	{
	$count_binomials2 ++;
	# print "$key\t" , $binomial_recorder{$key} , "\n";
	$count_binomials += $binomial_recorder{$key};
	}

print "\ncount_binomials:$count_binomials $count_binomials2\n\nfin.\n";
exit;







#####################################################################################################
#
#
#
#####################################################################################################



sub store_nodes
{

# nodes.dmp contains each taxon, described in a tree-like structure. 
# so each node has a name (eg polyphaga) and the higher group node to which it belongs (beetles). 

open(NODES , "nodes.dmp") || die "
ERROR. file nodes.dmp has not been found in current directory
instructions to get this file follow:
	please download (enter adress into a browser) the NCBI taxonomy database from:
	ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
	then unzip. it will create a new folder,
	enter the folder and copy the files 'names.dmp' and 'nodes.dmp' to the working directory.
	then your have the required file, and can try run this script again.
\n";

my $line_counter=0;
while (my $line = <NODES>)
	{
	if($line =~ /^(\d+)\t[\|]\t([^\t]+)\t[\|]\t([^\t]+)\t[\|]\t/)
		{
		my $tax_id = $1;my $parent_tax_id = $2;my $rank = $3;
		$rank_hash{$rank}++;
#		print "tax_id:$tax_id parent_tax_id:$parent_tax_id rank:$rank\n";

		if($tax_id == 48709){print "\n\n48709 stored\n"}

		my $current_rankcode="";
		 # if($rank eq "no rank"){$current_rankcode=""};
		 #if($rank eq "order"){$current_rankcode="";$current_order=$node->scientific_name();$current_suborder="";$current_superfamily="";$current_family="";$current_subfamily="";$current_tribe=""};
		 if($rank eq "suborder" ){$current_rankcode=1};
		 if($rank eq "infraorder"){$current_rankcode=2};
		 if($rank eq "series" ){$current_rankcode=2};
		 if($rank eq "superfamily" ){$current_rankcode=3};
		 if($rank eq "family"){$current_rankcode=4};
		 if($rank eq "subfamily"){$current_rankcode=5};
		 #if($rank eq "tribe" || $tribe_test=~ /^\w+ini$/){$current_rankcode=6;$current_tribe=$node->scientific_name();$current_tribe=~ s/ /_/g};
		 if($rank eq "tribe"){$current_rankcode=6};
		 # if($rank eq "subtribe" ){$current_rankcode=""};
		 if($rank eq "genus"){$current_rankcode="7"};
		 # if($rank eq "subgenus"){$current_rankcode=""};
		 # if($rank eq "species group"){$current_rankcode=""};
		 # if($rank eq "species"){$current_rankcode="";$in_species=1;}else{$in_species=0; # $current_species=$node->scientific_name()
				
		 # if($rank() eq "subspecies"){$current_rankcode="";$in_species=1};


		my $omit_current_taxon = 0;
		foreach my $omit_taxon(@omit_taxa)
			{
			if($omit_taxon == $tax_id){print "\nfound node ($omit_taxon) you have ordered to ignore\n";$omit_current_taxon = 1}
			}

		unless($omit_current_taxon == 1)
			{
			$nodes{$tax_id}{rank} = $rank;
			$nodes{$tax_id}{rank_code} = $current_rankcode;
			$nodes{$tax_id}{parent} = $parent_tax_id;
			$nodes{$parent_tax_id}{child_nodes} .= "\t" . $tax_id;
			}

		}else{
		print "line_counter:$line_counter line:$line";
		die "UNEXPECTED LINE:$line\nquitting\n";
		}
	$line_counter++;
	}

close(NODES);

my @ranks = keys %rank_hash;@ranks = sort @ranks;
print "ranks found in nodes.dmp:\n";
print LOG "ranks found in nodes.dmp:\n";

foreach my $rank(@ranks){print "$rank\t" , $rank_hash{$rank} , "\n";print LOG "$rank\t" , $rank_hash{$rank} , "\n"};

@all_nodes = keys %nodes;@all_nodes = sort @all_nodes;

# foreach $node(@all_nodes){print "$node\t" , $nodes{$node}{child_nodes} , "\n"};

}




#####################################################################################################
#
#
#
#####################################################################################################


sub parse_namesfile
{

# here just parse the scientific name of each node. ignore synonyms etc

open(NAMES , "names.dmp") || die "
ERROR. file names.dmp has not been found in current directory
instructions to get this file follow:
	please download (enter adress into a browser) the NCBI taxonomy database from:
	ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
	then unzip. it will create a new folder,
	enter the folder and copy the files 'names.dmp' and 'nodes.dmp' to the working directory.
	then your have the required file, and can try run this script again.

\n";

my $names_line_counter=0;
while (my $line = <NAMES>)
	{
	$names_line_counter++;#print "$names_line_counter\n";
# 24	|	Shewanella putrefaciens	|		|	scientific name	|

	if($line =~ /^(\d+)\t[\|]\t([^\t]+)\t[\|]\t([^\t]*)\t[\|]\tscientific name/)
		{
		my $tax_id = $1;my $name = $2;#my $rank = $3;
		# print "tax_id:$tax_id name:$name\n";

		# yet more strange chars, causing errors later:Apanteles_sp_BIOUG<CAN>_08BBHYM_1696_HQ552708


		# if you want to remove non-alphanumerical characters from assigned species names:
		$name =~ s/[\(\)\,\[\]\'\#\&\/\:\.\-\<\>]/ /g;
		$name =~ s/\s\s+/ /g;$name =~ s/\s+$//;$name =~ s/^\s+//;

		$nodes{$tax_id}{name} = $name;
		}else{
		if($line =~ /^(\d+).+scientific name/){die "UNEXPECTED LINE:\n$line\nquitting\n"}
		}

	}

close(NAMES);


}




#####################################################################################################
#
#
#
#####################################################################################################


sub traverse_nodes
{
my $current_node = $_[0];my $current_node_taxstring = $_[1];

# my $current_node = shift;

#if($current_node == 48709){
#print "\n\n$current_node $current_node_taxstring\n";
#}


my $child_nodes = $nodes{$current_node}{child_nodes};$child_nodes =~ s/^\t//;
my @child_nodes_array = split(/\t/, $child_nodes);

if($current_node == $starting_node)
	{
#	unshift @child_nodes_array , $starting_node;
	my $rank = $nodes{$starting_node}{rank};$rank =~ s/\s/_/g;$how_many{$rank}++;
	my $current_node_taxstring = substr $current_node_taxstring, 0,1;
	my $originalname = $nodes{$starting_node}{name};$originalname =~ s/[\s\t]/_/g;
	my $child_complete_lineage = $nodes{$starting_node}{complete_lineage} . "$rank:$originalname ";
	print OUT "$current_node_taxstring $originalname $starting_node $rank $child_complete_lineage\n";
	print OUT2 "$current_node_taxstring $originalname $starting_node $rank $child_complete_lineage\n";
	}




# print "\nchild_nodes:$child_nodes\n";die;

foreach my $child(@child_nodes_array)
	{
	unless($current_node == $child)# one the child nodes of the root node (1), is also 1 
	{
#	my $name = $nodes{$child}{name};
	my $rank = $nodes{$child}{rank};$rank =~ s/\s/_/g;$how_many{$rank}++;
	# my $current_rankcode = $nodes{$child}{rank_code};



	
	######################################################################################
        my $taxstring_and_name = get_taxstring( $child , $current_node_taxstring, $current_node);#
	######################################################################################
	
#	unless($parse_species_only ==1 )
#		{
	
		my @taxstring_array = split ( /__/, $taxstring_and_name );
		my $taxstring = $taxstring_array[0];
		my $name = $taxstring_array[1];
	
		my $name_copy = $name;$name_copy =~ s/\s/_/g;
	
		my $name_string = $nodes{$child}{name};$name_string =~ s/\s+/_/; ####################### sep2013


		my $child_complete_lineage = $nodes{$current_node}{complete_lineage} . "$rank:$name_string ";#prev assigned $name
		my $child_taxstring = $current_node_taxstring . $taxstring;
	
		$nodes{$child}{complete_lineage} = $child_complete_lineage;
	
		#print "complete lineage:$nodes{$child}{complete_lineage}\n";

	
		if(length($child_taxstring)>=30)
			{
		#	print LOG "NAME >=30chars:\nnode:$current_node child:$child name:$name child_taxstring:$child_taxstring lin:" , $nodes{$child}{complete_lineage} , "\n";
		
			}

#		}

#CA1D3D4C5Co6Co7pay Colymbetes_paykulli 183353 species suborder:Adephaga superfamily:Dytiscoidea family:Dytiscidae subfamily:Colymbetinae tribe:Colymbetini genus:Colymbetes species:Colymbetes_paykulli 

	my $originalname = $nodes{$child}{name};$originalname =~ s/[\s\t]/_/g;


my $name_assignment_to_taxnumber = "";

if($ignore_subspecies == 1)
	{
	if ( $rank eq "subspecies" || $rank eq "varietas"   || $nodes{$current_node}{complete_lineage} =~ / species:/)# 
		{
		my $parentoriginalname = $nodes{$current_node}{name};
		if($parentoriginalname =~ /^[A-Z][a-z]+\s[a-z]+$/)
			{
			$parentoriginalname =~ s/[\s\t]/_/g;
			$name_assignment_to_taxnumber = $parentoriginalname;
		#	print "node:$child named:$originalname rank:$rank appears to be subspecfic\n";
		#	print "\tassigning parent name instead:$parentoriginalname\n\n";
			}else{$name_assignment_to_taxnumber = $originalname}
		}else{
		$name_assignment_to_taxnumber = $originalname}

	}else{
	$name_assignment_to_taxnumber = $originalname
	}

print OUT "$child_taxstring $name_assignment_to_taxnumber $child $rank $child_complete_lineage\n";
print OUT2 "$child_taxstring $name_assignment_to_taxnumber $child $rank $child_complete_lineage\n";






#	unless($child == $starting_node)# only assign start node to a child for the purpose of printing taxonomic details of the root node to the output
#		{
		###########################################
		traverse_nodes($child , $child_taxstring);#
		###########################################
#		}
	}}
	
	
}






#####################################################################################################
#
#
#
#####################################################################################################


sub get_taxstring
{
my $current_node = $_[0];my $current_node_taxstring = $_[1];my $parent_node= $_[2];
my $name = $nodes{$current_node}{name};
my $rank = $nodes{$current_node}{rank};
my $rank_code = $nodes{$current_node}{rank_code};
my $newcode;
my $to_walk;

my $truncate_factor = 1;
if($rank eq "genus"){$truncate_factor =2};
if($rank eq "species" || $rank eq "subspecies"  || $rank eq "varietas" ){$truncate_factor =3};




# get rid of wierd chars

$name =~ s/\-//g;
$name =~ s/\'(\w+)\'/$1/g;
$name =~ s/\.//g;
$name =~ s/\://g;
$name =~ s/\// /g;
$name =~ s/[\(\)\,\[\]\'\#\&]/ /g;


my @words = split( / /, $name );

   if ( $rank eq "species" )
	{

if ( $name =~ /^[A-Z][a-z]+\s+[a-z]+$/){$binomial_recorder{$name}++}else{$nonbinomial_recorder{$name}  ++}

	$name=~ s/BMNH//;
        my @newwords;
        if ( $#words >= 2 ) #  THREE or more words
		{

            	foreach my $word (@words) 
			{
                	unless($word =~ /^sp$|^aff$|^nr$|^cf$|^spp$|^$/ ){push( @newwords, $word )};
            		}

		$to_walk = $newwords[$#newwords];if($to_walk =~ /\s/){die "\nerror 391\n"}
		if($to_walk =~ /\d/ && length($to_walk)<=10){$truncate_factor = length($to_walk)};
		$newcode = walk_name(	$to_walk , $current_node_taxstring , $truncate_factor , $current_node);

        	}else 
		{

	        if ( $#words == 1 ) #  TWO words
			{
			
			$to_walk = $words[1];
        		$newcode = walk_name(	$to_walk , $current_node_taxstring , $truncate_factor , $current_node);
			#if($to_walk =~ /^sp$|^aff$|^nr$|^cf$|^spp$/ ){$unidentified_sp_hash{ $words[0] . "__" . $words[1]}++ }
			
			}else{
			
			$to_walk = $words[0];
			$newcode = walk_name(	$to_walk , $current_node_taxstring , $truncate_factor , $current_node);
			
			}
			
        	}

	}
	elsif ( $rank eq "subspecies" || $rank eq "varietas"   || $nodes{$parent_node}{complete_lineage} =~ / species:/)# 
	{
 
#	unless($ignore_subspecies==1)
#		{

		$to_walk = $words[$#words];
		$newcode = walk_name(	$to_walk , $current_node_taxstring , $truncate_factor , $current_node , 1);
#		}

	}else{	#if not species or below

	# if you want to shorten codes, could tell it to ignore atypical ranks here, such as subgenus, species complex, no rank, etc

#	 print "passing n:$name  c:$current_node_taxstring to subwalkname\n";



	unless($parse_species_only == 1)
		{
	if($ignore_atypical_ranks_to_shorten_codes==0)
		{

		if($name =~ /\s\S/){
			print LOG "\n>1 word for non-species rank -> ($name)";$name =~ s/.+\s(\S+)/$1/;print LOG "\n\tusing last word as name ($name)\n"}

		$to_walk = $name;
	        $newcode = walk_name($to_walk , $current_node_taxstring , $truncate_factor , $current_node);
		}else{

	#	if( $rank eq "varietas" || $rank eq "subspecies" || $rank eq "species" || $rank eq "genus" || $rank eq "tribe" || $rank eq "subfamily" || $rank eq "family" || $rank eq "superfamily" || $rank eq "suborder" || $rank eq "order" )
	#		{
			if($name =~ /\s\S/){print LOG "\n>1 word for non-species rank -> ($name)";$name =~ s/.+\s(\S+)/$1/;print LOG "\n\tusing last word as name ($name)\n"}
			$to_walk = $name;
        		$newcode = walk_name($to_walk , $current_node_taxstring , $truncate_factor , $current_node);
	#		}

		}

		}
	}

	my $newcode_plus_name = $newcode . "__" . $to_walk;
#	print "\t$rank $current_node_taxstring newcode$newcode\n";
	return ($newcode_plus_name);

}




#####################################################################################################
#
#
#
#####################################################################################################


sub walk_name
{
my $to_walk = $_[0];my $current_node_taxstring = $_[1];my $truncate_factor= $_[2];
my $current_node = $_[3];my $subsp_test = $_[4];
my $name = $nodes{$current_node}{name};
my $current_rankcode = $nodes{$current_node}{rank_code};

# input is taxnomic name at current node

my $current_index = 0;
my $test_this_string = "";
	

if($ignore_subspecies==1 && $subsp_test == 1)
	{

	}else{

	$test_this_string = $current_node_taxstring . substr($to_walk , 0,  ($truncate_factor+$current_index));

	while(exists($all_taxstrings{$test_this_string}))
		{
		$current_index++;
		if($current_index<=100)
			{
			if($duplicate_taxonstring_numbers == 1 && $truncate_factor == 1)
				{
				$test_this_string = $current_node_taxstring . $current_index . substr($to_walk , 0,  $truncate_factor);
				}else{
				$test_this_string = $current_node_taxstring . substr($to_walk , 0,  ($truncate_factor+$current_index));
				}
			}else{
			if($test_this_string =~ /^(.+[_])(\d+)$/){$test_this_string = $1 . ($2+1)}else{$test_this_string = $test_this_string . "_1" }
			}
		};


	$all_taxstrings{$test_this_string} = 1;

	# my $to_return = substr($to_walk , 0,  ($truncate_factor+$current_index));
	# my $to_return = $test_this_string; NO! this is the total string, we need to return just the new substring

#	print "test_this_string:$test_this_string\n";
	$test_this_string =~ s/$current_node_taxstring(.+)$/$1/;


	if($rank_numbers_in_taxcodes == 1)
		{$test_this_string .= $current_rankcode}

	}




return($test_this_string);

}


#####################################################################################################
#
#
#
#####################################################################################################


sub read_command_arguments
{

my $arguments = shift;
$arguments =~ s/\n//;$arguments =~ s/\r//;

print "
\n
\n";



# $starting_node = $ARGV[0];
if($arguments =~ /-starting_node\s+(\S+)/)
	{
	$starting_node = $1;
	}else{
	die "\nerror. please give starting node in command, as so: -starting_node [ncbi_taxon_number]\n\n";
	};

if($arguments =~ /-ignore_subspecies/)
	{
	$ignore_subspecies = 1;
	}else{
	$ignore_subspecies = 0;
	}
# $ignore_subspecies			= 1;	# 0=build codes beyond species level (e.g. include subspecies, varietas, and unnamed bacteria ranks)
						# if 1, these extra ranks are ignored, and species strings are assigned to them 
						# this applies to both tobycode and full species name, assigned to given ncbi tax number


if($arguments =~ /-parse_species_only/)
	{
	$parse_species_only			= 1;	# default = 0
	}else{
	$parse_species_only			= 0;	# default = 0
	};
# $parse_species_only			= 1;	# default = 0

# -starting_node -ignore_subspecies -parse_species_only

print "

user has specified:
	starting_node:$starting_node
	ignore_subspecies:$ignore_subspecies
	parse_species_only:$parse_species_only

";


};

#####################################################################################################
#
#
#
#####################################################################################################








