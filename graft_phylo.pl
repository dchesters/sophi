#!/usr/bin/perl




###########################################################################################################################################
#	
#	
#	
#	
#	
#	relational_constraints.pl / graft_phylo.pl
#		  
#    	Copyright (C) 2015-2022  Douglas Chesters
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
##########################################################################################################################################
#
#	
#	
#	
#	Notes as of 2019-AUG-07: current version of the whole process is in a relativly good condition. Will finalize as is.
#		This script is one of a 4 step process:
#		1) read taxonomic database and produce a processed file ready for input.
#		2) read backbone phylogeny infer taxon constraints and write in fasttree format and as characters that can be read into raxml
#		3) herein, read backbone phylogeny and infer relational constraints, 
#			write as fasttree format, and as newick string as a regular raxml constraint tree.
#	
#	
#	
#	
#	
#	
#	
#	
#	CHANGE LOG
#
#	21 Jan 2015: 	Version 1.001
#	22 oct 2015: 	Modified some bacterial labels when reading ncbi taxonomy db:/^Candidatus\s(\w+)/$1/
#			these are often not used, at least in the tree  
#			if cant find taxonomy for binomials, tries genus name
#	01 dec 2015: 	Runs without outgroup specified
#			Finds taxonomies for genus name when lineage for species is not found
#	28 Jan 2018: 	Implemented another algorithm,
#			an advance over the 'basic constraint method'
#			makes newick style constraint tree, will usually omit non-overlapping taxa,
#			maintains backbone topology, inferrs higher ranks for backbone terminals to the extent possible
#	22 Mar 2018: 	More helpful error message when user uses out of date ncbi taxonomy database
#			got it working on backbone trees with just genus names as labels
#	01 Nov 2018: 	User specifies ranks at which constraints are made, in command.
#			avoids 'bug', early version was hacked so order constraints not made,
#			for reasons that might make sense at the time,
#			but which meant inadequate constraining since.
#	21 Dec 2018: 	Prints fasttree format constraints of same set printed to newick.
#			although exsiting scripts available for fasttree format relational constraints,
#			this script might work better in cases of low overlap and missing taxonomies.
#			renamed, backbone_constraints_newick.pl to relational_constraints.pl
#	17 Jul 2019: 	Taxonomic database parsing moved to seperate script (taxon_table.pl),
#			allows user to check taxonomies and whether anything of particular interest is satisfactory,
#			and allows manual modification.
# 			Incorporated ITIS taxonomic database.
#			Discontinued user setting of constraint ranks, it was interfeering with reading of lineages.
#			would need to be done at a later point.
#	18 Jul 2019:	Error message if no taxon overlap between tax table and phylogeny, indicating check tax number input to previous script 
#			*** Fixed Newick based relational functions. now automatically prunes constraint tree to only barcode taxa, 
#			vastly easier for user in cases of low overlap between backbone and barcodes. 
#			Code is integrated into fasttree sub for equivelency, makes exactly same set of constraints.
#	06 Oct 2019: 	Discontinued selected definition of ranks (option -constrain_ranks, and even default rank list), 
#			was making strangly incomplete lineages for many terminals.
#	05 Mar 2020: 	Error message if users backbone tree has unparsable terminal labels
#	07 Mar 2021: 	Incorrect output file was being used by some users, so correct one is printed to screen as final message			
#	10 Mar 2021: 	Prints list of fasta IDs which are constrained (was an earlier list but might not be exact)
#	06 Jun 2021: 	Added function to plot image depicting backbone tree and branches applied, with indicators of counts for correponding barcode
#	20 Sep 2021: 	New outfile (lineage_assignments_to_barcodes) as need to know in later processes exactly how this script assigned taxa.
#	11 Oct 2021:	Minor change to labelling in r image commands, now gives tax+counts 
#	12 Oct 2021:	Code for making R plot of constraint process now as user option (-plot_constraints)
#				plotting requires removal of branchlengths (due to quick writing), which is done if this option is used
#				(was uncertain how other functions would be affected by branchlength removal, so wont do that if this option not specified)
#	04 Jul 2022:	A few improvements of R visualization of constraint process, such as option to print terminal labels on backbone	
#	28 Jul 2022:	Major new function; grafting source trees (or components thereof) to equivelent terminals of single backbone.
#	06 Aug 2022: 	Improved function for assigning taxa to node; doesnt break so often where there is less higher taxon retrieval,
#			which happens a lot on morphological source trees for which genus names are not in ncbi tax db.
#			Reduced number of taxon/clades stored for each source tree in grafting algorithm
#			(dont bother storing one encountered if already stored a larger one for that taxon)
#			Additional check for grafting cutoff.
#	08 Aug 2022:	Minor fix, dont open file grafting log if not using that function.
#	20 Aug 2022:	Minor fix of grafting alrorithm.
#	22 Aug 2022:	Option to retain graft node labels through multiple graft iterations (to be used for plotting). 
#	23 Aug 2022:	Bugfix, unneccesary pair of parentheses added when grafting.
#	26 Aug 2022:	Minor changes.
#	16 Sep 2022:	Option to exclude a terminal from inference and replacement with broader sample,
#			Probably want to do this if single outgroup is being used.
#	16 Sep 2022:	As grafting (implemented first in relational_constraints.pl) has become a key component in expanding backbone topologies,
#			will fork this script to seperate grafting from downstream step of consolidating a backbone with species-rich DNA matrices.
#			Forked script will be named graft_phylo.pl
#	
#	
#	
#	
#	
#	
#	
#	
##########################################################################################################################################

#  OLD notes (may no longer apply)




# as opposed to selecting a shared taxon, there would often be a few ranks above that that node would represent,
# need to draw these on, then select a criteria for assigning, mybe there would be ambiguity.
# 
# 
# second type of constraint made, monophyly of the assigned name, 
# at first instance root->tip, if contains no other taxa, and >1 of the taxa of rank below.
# 
# 
# 
# 
# SUBS:
# read_command_arguments, store_nodes, parse_namesfile, traverse_nodes, read_fas, record_tree2
# traverse_backbone_tree_and_infer_taxonomic_node_labels {get_terminals_from_this_node2, get_shared_taxa, get_all_taxa}
# traverse_backbone_tree_and_collapse
# 
# 
# 
# 
# 
# 
# 
# sub store_nodes reads file nodes.dmp, which contains ncbi taxonomic heirachy, 
# 	with nodes represrted by unique numbers
# sub parse_namesfile parses names.dmp for 'scientific name' of each node
# sub traverse_nodes, recurses through all nodes of the taxonomic heirachy, from user specified position
# 	constructs lineage inforamtion for each species
# sub read_fas, reads users fasta file, stores all taxa included,
# 	retrives lineage for each, and stores all higher taxa.
# sub record_tree2, read user tree in newick format,
# 	some basic processing if required, stores structure of tree,
# 	reads tip taxa, again rretreives their higher taxa and records these.
# sub backbone_constraints2
# 	for each terminal of teh backbone tree,
# 	retrives lineage, then proceeds through all ranks present, starting from genus, in a root-ward direction  
# 	higher taxon names are assigned to the terminal, until the point is readch that the higher taxon name 
# 	is possessed also by a terminal anywhere else on the tree. 
# 	note, in the unlikely event of complete monophylys on the backbone tree,
# 	it would be neccessary only to consult the taxa of the terminals sister to the termnial being assigned.
# the sub prints a newick string just of these assignments, for information.
# next , for the taxon at each terminal, fasta IDs belonging to that terminal are retrived ,
# and , where >1, are joined in comma-sepeatred parenthseis,
# and inserted into the newick string , to replace its encompasing taxon name
# finally, any terminals from taxa which are absent in the sequence data,
# are pruned from the constraint tree.
#
#
#  
# 
# specifications for phylogeny: directional, preferably unrooted. since raxml constraint tree input requires unrooted.
# 				absolutly no duplicate tip IDs. perferably binomial Genus_species
# 
# 
# 
# 
# fork of read_deep_level_constraints.pl
# doesnt set constrAINTS that require floating taxa (see fasttree constraints),
# so can be set in strict newick format
# 
# 
# 
# 
# 
# 
# 
# 
# 
#################################################################################################################################


my $arg_string  = join ' ', @ARGV;



# additional user options:

$verbose 		= 1; # 0 simple screen output, 1=verbose, 2=ulta-verbose

#$tree_ID_format = 2;# i think 1=Order_Genus, 2= Order_Genus_species


# $process_backbone_tree_terminal_IDs = 2;
# if == 1, then process ids in backbone tree, e.g. Coleoptera_Apatides_fortis to Apatides_fortis
# if == 2, then binomial to genus

$fasttree_format_constraints_file 	= "ftcons"; # output.
$fasttree_constraints_tabdelim 		= "relational_constraints_tabdelim";


$process_new_consraint = 0; # default 0, 
				# somtimes there will be taxa inferred from the input backbone tree
				# which are not in the subsequent dataset,
				# usually rare, but if == 1 , script will try prune them from the 
				# subsequent constraint tree


$ugly_hack = 1;

$unconstrain_broadly_labelled_data = 1;
# much data only ident for order level, usually ignored, although many slip through, this is an extra chaeck
# it tests if the species has a family designation, if not, then not constrained.


			# currently only one option implemented
$newick_reformat = 1;# 1 = backbone ID are binomial, convert to genus name only


			# for R plot of rel constraint process
$barcode_box_taxname_cutoff 		= 5;# cutoff for printing tax name for a barcode box.
$terminal_labels_on_R_visualization 	= 1;
$r_visualization_barcode_box_tax_cex 	= 0.5;

$grafting_cutoff 	= 5; # appies only during grafting option, if subtree of matching taxon is below this, wont be grafted


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
$outgroup = "NULL";
$outfile_prefix;



#####################################
read_command_arguments($arg_string);#
#####################################

$output_filename = "$fas_file.$treefile";


###############################################################################################################################################


# sub store_tax_heirarchy stores the following items:
# 	ncbi_nodes{}{rank_code}; $ncbi_nodes{}{rank}; $ncbi_nodes{}{parent}; 
#	$ncbi_nodes{}{child_nodes}; $ncbi_nodes{}{name}

#######################
store_tax_heirarchy();#
#######################


#################################################

# traverse nodes starts at first parent node in the input file, which should be the root.
# the input file was made with a root node user specified.


# sub traverse_nodes uses:
#  ncbi_nodes{}{child_nodes}; $ncbi_nodes{}{rank};$ncbi_nodes{}{name};$ncbi_nodes{}{rank_code};
# and writes:
#  complete_lineage_for_this_species ; ncbi_taxnumber_for_taxname	

$starting_name = $ncbi_nodes{$starting_node}{name};$starting_name_ff = $starting_name;$starting_name_ff =~ s/\s+/_/g;
print "name for the starting node ($starting_node)\nwhich has been specified ($starting_name)\n";
$starting_name =~ s/^(\w).+$/$1/;

#################################################
traverse_nodes($starting_node , $starting_name);#
#################################################


my @ranks_array = keys %how_many;
print "heirarchy traversed:
	nodes_traversed:$nodes_traversed
	ranks encountered:$#ranks_array
";



###############################################################################################################################################


open(BCN_LOG, ">backbone_constraints_newick_LOG") || die "";
close BCN_LOG;


my %terminals = ();
my $root_identity = "";
my $date = localtime time;

%count_tax_in_new_data;
%query_IDs;	

 # the fasta file contains names of new entries to be placed onto backbone tree. 
 # current limitation that fasta IDs need to be in the format >Genus_species (e.g. Aenictus_aratus)
 # 	and have no duplicate IDs (ie filtered at the species level)
 # fasta IDs are stored in the hash %query_IDs.

print "\ncalling sub read_fas\n";
open(BCN_LOG, ">>backbone_constraints_newick_LOG") || die "";

############
read_fas();#
############

close BCN_LOG;

####################################################################################################################################

# global
$current_sourcetree;

if($graft_source_trees_list =~ /\w/)
	{
	print "\ncalling sub read_source_trees\n";
	read_source_trees();
	};

####################################################################################################################################


%nodes;		# user phylogeny stored in this hash
$root_node;	# this global will keep the value of the last in the loop, which will be the root
%count_taxa_represented_in_terminals;	# will use this to determine whether taxa are monophyletic.
$print_internal_nodes=0;

if($retain_backbone_node_labels == 1)
	{

	#########################
	record_tree_graft($treefile);#
	#########################
	
	}else{

# read backbone tree into hash. node identities for hash key, 
# then parent/child identities as hash entry
# gets complete lineage for each tip, and for tax at all ranks, 
# counts number of tips for each into object count_taxa_represented_in_terminals
# optional processing of tip labels. otherwise, this sub just does bog-standard reading of tree

# as a raxml constraint tree will be made, here the script will try and unroot if a rooted tree is input

#########################
record_tree2($treefile);# 	$treefile is the backbone tree
#########################

	};

my @array_taxa_found_on_backbone = keys %count_taxa_represented_in_terminals;

print "
your backbone tree has been read. 
	TIPS:$count_the_terminal_nodes 
	ALL NODES:$print_internal_nodes.
	total taxa found on backbone: $#array_taxa_found_on_backbone
\n";




####################################################################################################################################

# jan 2018 , another constraint method: Bottom-Up

open(BCN_LOG, ">>backbone_constraints_newick_LOG") || die "";

#########################
backbone_constraints2();#
#########################

close BCN_LOG;

# Sub contains everything, output file now printed
# Important output now printed to file:$outfile_prefix.less_basic_constraint_tree
# Also informative is $outfile_prefix.less_basic_constraint_tree.internal_tax_labeled

if($graft_source_trees_list =~ /\w/)
	{
	print "\n\nUser opted for grafting algorithm, done. Made $grafts_made seperate subtree grafts, totalling $grafted_count terminals\n";
	print "\toutput files: backbone_with_grafted_sourcetrees, backbone_with_grafted_sourcetrees.internal_labelled\n\n";
	exit;
	};

#############################################################################################

# 2021-06: For visualization. This function is transferred from process_newick. 
#		It makes counts of terminals from each node, used in calculating y positions.

print "Traverse backbone topology 1, counting terminals from each node, to be used later in visualization\n";

$traversal_start_node = $root_node; 

#######################################################################
assign_names_to_internal_nodes($traversal_start_node , "New_Root");#
#######################################################################

print "Fininshed traverse backbone topology 1\n";

#############################################################################################


get_sum_branchlength($root_node , 0, 1);
print "maximum_x will be:$maximum_x\n";
$plot_y_min = 1000000;


# AUG 2019:
#	 try some new code, it is awkward to prune extranous stuff off the constraint tree later, particularly if low overlap backbone was used.
#	see if that can be done here.
#	conduct from within the fasttree section, this is already retreiving barcode taxa terminals for each backbone node




print "\nprinting fasttree format relational constraints\n";
$ft_relational_constraints_subcalls =0;
$pruned_constraint_newick = "($root_node)";# $pruned_constraint_newick is the main output constraint tree

# **** NOTE **** this function produces the main output file, newick constraint tree with extranous stuff ommited



#####################################################################
# also plotting commands
$xlim 				= "c(0,3)"; $ylim 	= "c(0.4,1.6)";
$image_settings_RCT 		= "tiff\(filename = \"ITOL3_overlap.tiff\", width = 8000 , height = 9000, units=\"px\", res=600, compression = \"jpeg\"\)\n";

if($plot_constraints == 1)
	{
	open(OUT_RECTANGL , ">r_commands_filename_RECTANGL") || die "\nerror 241. cant open output file for writing r script ($r_commands_filename)\n";
	print OUT_RECTANGL "library(plotrix)\n$image_settings_RCT";
 	print OUT_RECTANGL "colorfunc = colorRamp(c(\"blue\" , \"gray\", \"gray\", \"gray\", \"gray\", \"gray\", \"gray\", \"gray\", \"gray\", \"gray\", \"gray\"))\n";
	print OUT_RECTANGL "plot(0, type = \"n\", xlim = $xlim, ylim = $ylim,";
	print OUT_RECTANGL " xlab = \"\",ylab = \"\",xaxt = \"n\", yaxt = \"n\", bty = \"n\")\n";
	};

#####################################################################




###############################################################
print_fasttree_format_relational_constraints($root_node , 0, 1);# # $root_node takes value of last node read of the backbone tree, the root
###############################################################

 # now stored in %relational_contraints


print "left sub\n
sucessful_constraint_insertions_A:$sucessful_constraint_insertions_A
sucessful_constraint_insertions_B:$sucessful_constraint_insertions_B

";


@terminal_ys = keys %corresponding_barcodes; @terminal_ys = sort { $a <=> $b } @terminal_ys;
my $total_barcodes_plotted;
foreach my $terminal_y(@terminal_ys)
	{
	my $count_barcodes_this_terminal = $corresponding_barcodes{$terminal_y};$total_barcodes_plotted += $count_barcodes_this_terminal;
	if($count_barcodes_this_terminal >= $maxbarcodes){$maxbarcodes = $count_barcodes_this_terminal};
	}

my $current_y1 = 0;
foreach my $terminal_y(@terminal_ys) # these are terminals
	{
	my $count_barcodes_this_terminal = $corresponding_barcodes{$terminal_y}; # and hash contents are barcodes assigned to backbone terminal
	my $current_y2 = $current_y1 + $count_barcodes_this_terminal;
	my $taxY = $tax_of_Y{$terminal_y};


	my $y1_adjust = $current_y1 / $total_barcodes_plotted; my $y2_adjust = $current_y2 / $total_barcodes_plotted; 
	$y1_adjust += 0.5;$y2_adjust += 0.5;
#	my $y_diff = $y2_adjust - $y1_adjust;
	my $y_mean = ($y2_adjust + $y1_adjust) / 2;
	my $branchthickness = 0.5;my $barcode_color = "gray";

	my $R_command =	"scaled_val <- ( $count_barcodes_this_terminal - 1 ) / ( $maxbarcodes - 1 ); " . 
			"colors2<- rgb(colorfunc(scaled_val) , maxColorValue=255)\n";
	$draw_tree_R_commands_RECTANGL .= $R_command;



	# Text, taxonomic name
	if($count_barcodes_this_terminal >= $barcode_box_taxname_cutoff)  
		{
		if($taxY =~ /\w/)
			{
			my $label6 = $taxY . " " . $count_barcodes_this_terminal;
			my $R_command = "text(2.7,$y_mean,adj=c(0, 0.5),labels=\"$label6\",cex = $r_visualization_barcode_box_tax_cex)\n";$draw_tree_R_commands_RECTANGL .= $R_command;
			};
		$branchthickness = 1;
		};

	# LINE
	if($count_barcodes_this_terminal >= 20)
		{
		$barcode_color = "black";my $R_command = "segments( 1, $terminal_y , 2, $y_mean , col = \"black\", lwd = 0.25)\n";	$draw_tree_R_commands_RECTANGL .= $R_command;
		};


	# colored BOX
	# lty=NULL, , xpd=FALSE
	my $R_command =  "rect(2, $y1_adjust, 2.5, $y2_adjust, col=colors2, border=NA, lwd=0.5)\n"; # border=\"black\"
	$draw_tree_R_commands_RECTANGL .= $R_command;

#	my $R_command =  "segments(2, $y_mean, 2.5, $y2_adjust, col=\"$barcode_color\", lwd=$branchthickness)\n";
#	$draw_tree_R_commands_RECTANGL .= $R_command;
#	my $R_command =  "segments(2, $y_mean, 2.5, $y1_adjust, col=\"$barcode_color\", lwd=$branchthickness)\n";
#	$draw_tree_R_commands_RECTANGL .= $R_command;



	$current_y1 = $current_y2;
	}


if($plot_constraints == 1)
	{
	print OUT_RECTANGL $draw_tree_R_commands_RECTANGL;# segments, points, text
	print OUT_RECTANGL "dev\.off\(\)\n";
	close OUT_RECTANGL;
	print "run this command to visualize: R < r_commands_filename_RECTANGL --vanilla --slave\n";
	};



# pruned_constraint_newick:$pruned_constraint_newick

# not usually needed:
# open(OUT_CPP, ">constraint_phylogeny_pruned.nwk") || die "\nerorr 404\n";
# print OUT_CPP "$pruned_constraint_newick;\n";
# close OUT_CPP;

$pruned_constraint_newick =~ s/^\(//;$pruned_constraint_newick =~ s/\)$//;
open(OUT_CPP2, ">constraint_phylogeny_pruned2.nwk") || die "\nerorr 404\n";
print OUT_CPP2 "$pruned_constraint_newick;\n"; # $pruned_constraint_newick is the main output constraint tree
close OUT_CPP2;

my @anotherlist = keys %another_list_of_RelConstrained_barcodeIDs;@anotherlist = sort @anotherlist;
open(OUT_C4, ">list_constrained_IDs") || die "\nerorr 398\n";
foreach my $ID(@anotherlist)
	{
	print OUT_C4 "$ID\n";
	}
close OUT_C4;

print "

 *** READ THIS! *** The important output file that you probably need is now printed, named:constraint_phylogeny_pruned2.nwk

 Note, under normal circumstances you DONT NEED to prune anything (input or output), script should account for overlap on its own.
 A list of sequence IDs that could be assigned to constraint tree in file: list_constrained_IDs. 
 If the following steps are slow, you can probably just skip (CTRL-C).
";



# relational constraint for each relevent backbone node, now in this object:
my @rel_constraints = keys %relational_contraints;@rel_constraints = sort @rel_constraints;

# something wierd happening, bunce of hash keys of fasta ID appear as 1.
# cant fingure out why, jsut read in again:
open(IN, $fas_file) || die "\nerror ... fasta file you have given ($fas_file) cannot be opened. \n";
while (my $line = <IN>)
	{
	$line =~ s/\n//;$line =~ s/\r//;
	if($line =~ /^>(.+)/)
		{
		my $id = $1; if($id == 1){die "\nWTF\n"};
		$fasta_IDs_yet_again{$id} = 1;
		}
	}
close IN;


my @the_fasta_IDs = keys %fasta_IDs_yet_again; # bugfix 201908
@the_fasta_IDs = sort @the_fasta_IDs;

	foreach my $id5(@the_fasta_IDs)
		{
	#	print "id5:$id5 ";
		};
# die "";

my $printing_limit = 0;
for my $i(0 .. $#rel_constraints)
	{
	my $rel_constraint_ID = $rel_constraints[$i];
	if($printing_limit <= 6){print "FT relational constraint $i of $#rel_constraints ($rel_constraint_ID)\n"};
	if($printing_limit =~ /00$/){print "FT relational constraint $i of $#rel_constraints ($rel_constraint_ID)\n"};

	$printing_limit++;

	my $rel_constraint = $relational_contraints{$rel_constraint_ID};
#	print "\ni:$i rel_constraint_ID:$rel_constraint_ID contraint:$rel_constraint\n";
	my @split_csontr = split /RelationalConstraint\d/ , $rel_constraint;
	# RelationalConstraint0

	foreach my $id(@the_fasta_IDs)
		{
		unless($id == 1) # im stumped why bunch of 1's appeared and persisting as hash keys,
		{		# edit 201907, seems because keys command missing when defining array

	#	print "id:$id\n";
		my $state = "-";
		foreach my $j(0 .. $#split_csontr) 
			{
			my $contraint7 = $split_csontr[$j];
			if($contraint7 =~ /\t$id\t/){$state = $j};
			};
		if($state =~ /\d/){$count_states_each_member{$id}++};
		$relational_constraints_per_fasta_ID{$id} .= $state;
		$relational_constraints_per_fasta_ID2{$id} .= "$state\t";
		if($id == 1){die "\nHUH!\n"}
		};
		};

	};

##################

print "\nlength array fasta_IDs: " , scalar @the_fasta_IDs , "\n";

open(OUT_REL_CONS, ">relational_constraints.ft_format") || die "";
open(OUT_REL_CONS2, ">relational_constraints.ft_format.tdt") || die "";
$printystring = "member\t";
for my $i(0 .. $#rel_constraints)
	{
	my $rel_constraint_ID = $rel_constraints[$i];$printystring .= "$rel_constraint_ID\t";
	};
print OUT_REL_CONS2 "$printystring\n";

foreach my $id(@the_fasta_IDs)
	{
	unless($id == 1) # bug fixed , this should not be neccessary now 
	{
	my $current_line = $relational_constraints_per_fasta_ID{$id};
	my $current_line2 = $relational_constraints_per_fasta_ID2{$id};
	my $fam = "NA";
	if($family_of_fasta_ID{$id} =~ /\w/){$fam =  $family_of_fasta_ID{$id}};

	print OUT_REL_CONS "$id\t$current_line\n";
	print OUT_REL_CONS2 "$fam $id\t$current_line2\n";
	}};

close OUT_REL_CONS;
close OUT_REL_CONS2;

foreach my $id(@the_fasta_IDs)
	{
	if($count_states_each_member{$id} =~ /\d/){$members_with_ft_rel_constraint++};
	};

my $count_ft_rel_constraints = scalar @rel_constraints;
my $count_fasta_ids = scalar @the_fasta_IDs;

print "
for fasttree relational constraints,
 number of constraints:$count_ft_rel_constraints
 number of members in fasta file:$count_fasta_ids
 of which, number with at least one constraint:$members_with_ft_rel_constraint
";

my @backbone_terminals_for_which_taxonomy_are_not_found = keys %backbone_terminals_for_which_taxonomy_not_found;

if(scalar @backbone_terminals_for_which_taxonomy_are_not_found >= 1)
	{
print "
\n\tWARNING   .....   taxonomic information not found for the following backbone terminals:\n\t\t\t@backbone_terminals_for_which_taxonomy_are_not_found
	recommend you prune these and rerun (behaviour uncertain if these are left in)
";
	};


print "\n\noutput file you probably need is: constraint_phylogeny_pruned2.nwk
	list of taxa in this constraint tree: list_constrained_IDs
";

print "\n\nFIN.\n\n";
exit;



####################################################################################################################################







%all_taxa_on_backbone_tree;
%monophylys;

	# traverse nodes of the user genome tree, 
	# for use in constraining treesearch in the barcode data


$new_newick_string = "($root_node)";
$new_newick_string2 = "($root_node)";
$newick_print_string = "($root_node)";



print "\n
calling sub traverse_backbone_tree_and_infer_taxonomic_node_labels, 
which will traverse nodes of the backbone tree, 
infer taxonomic names for each node
";

# this code maintained in taxon_constraints

#########################################################################
traverse_backbone_tree_and_infer_taxonomic_node_labels($root_node , 0);#
#########################################################################

my @alltax9 = keys %all_taxa_on_backbone_tree;@alltax9 = sort @alltax9;

open(NEWICK1 , ">$outfile_prefix.newick_string1") || die "\nerror 193.\n";
print NEWICK1 "$new_newick_string\n";
close NEWICK1;
open(NEWICK2 , ">$outfile_prefix.newick_string2") || die "\nerror 193.\n";#
print NEWICK2 "$new_newick_string2\n";
close NEWICK2;


my @mps = keys %monophylys;@mps = sort @mps;
foreach my $mp(@alltax9)
	{
	if(exists($monophylys{$mp}))
		{
		print "$mp IS monophyletic\n";
		}else{
		#print "$mp NOT monophyletic\n";
		
		}
}
print "\n$#mps name monophylys found on backbone tree, these will be set for next tree search if relevent 
, in addition to the other things\n";


# taxonomic node labels assigned to backbone tree.
# now traverse tree again and decide at which point to stop (collpase all decendents)
# this must be done where nodes are reached at which not all taxon represented in the fasta file 
# are included as decendents.
# also might need to collapse at where polyphyletic taxa are reached.



print "\ncalling sub traverse_backbone_tree_and_collapse
which will do something.
";


$yet_another_newick_string = "($root_node)";


#########################################################################
traverse_backbone_tree_and_collapse($root_node);#
#########################################################################


#while($yet_another_newick_string =~ s/\,*\(\)//){};
#while($yet_another_newick_string =~ s/\(\,\(/"(("/e){};
while($yet_another_newick_string =~ s/(\,)\,+/$1/){};
while($yet_another_newick_string =~ s/\(\,\)//){};
while($yet_another_newick_string =~ s/\(\,\(/"(("/e){};
while($yet_another_newick_string =~ s/\(\,([A-Z])/($1/){};
while($yet_another_newick_string =~ s/\,\)/)/){};


	# e.g. $yet_another_newick_string = ((Lepismatidae-Nicoletiidae)Zygentoma);
open(NEWICK2000 , ">$outfile_prefix.newick_string2000") || die "\nerror 193.\n";#
print NEWICK2000 "$yet_another_newick_string;\n";
close NEWICK2000;



# this object contains all names of taxa written onto the constraint tree:
# constraints_written


my @constraints_writtenn = keys  %constraints_written;@constraints_writtenn = sort @constraints_writtenn;
foreach my $constraint(@constraints_writtenn)
	{
#	print "constraint:$constraint\n";
	my @split = split /\-/ , $constraint;
	foreach my $item(@split)
		{
	#	print "\t$item\n";
		$all_constraint_taxa{$item} =1;
		};
	};





#my @ft_constraints_again = keys  %ft_constraints;
#@ft_constraints_again = sort @ft_constraints_again;
#open(CONS_TABLE , ">$fasttree_constraints_tabdelim") || die "\nerror 258\n";
#foreach my $constraint(@ft_constraints_again)
#	{
#	#print "constraint:($constraint)\n";
#	$constraint =~ s/\n/\t/g;#
#	$constraintno++;
#	print CONS_TABLE "constraint:$constraintno\t$constraint\n";
#	}
#close CONS_TABLE;




my @new_terminals = keys %query_IDs;@new_terminals = sort @new_terminals;

# my @all_taxa = keys %all_taxa_on_backbone_tree;@all_taxa = sort @all_taxa;


# open(OUT5 , ">$fasttree_format_constraints_file") || die "\nerrro\n";
# print OUT5 " " , scalar @new_terminals , " " , scalar @ft_constraints_again + scalar @mps + 1, "\n";

# fasttree constraint file.
# 5 2
#A        00
#B        0-
#C        10
#D        11
#E        -1

# comorpha)Brachycera)Diptera)Diptera-Siphonaptera-Mecoptera)Diptera-Siphonaptera-Mecoptera-Amphiesmenoptera,((S


my $coutn = 0;

open(OUT_LIST, ">$outfile_prefix.list_of_constrained_members") || die "\nerror 410\n";

foreach my $new(@new_terminals)
	{

	$coutn++;
	if($coutn=~ /0000$/)
		{
		print "$coutn of $#new_terminals, ID:$new\n";
		};
	my $genus = $new;$genus =~ s/^([A-Z][a-z]+)_.+/$1/;
	
	my $lineage = "NA";

	if(exists($complete_lineage_for_this_species{$new}))
		{
		# complete_lineage_for_this_species, not just species, 
		# but has lineage for every taxon in ncbi, only mod is space->underscore
		$lineage = $complete_lineage_for_this_species{$new};

		}elsif(exists($complete_lineage_for_this_species{$genus}))
			{
			$lineage = $complete_lineage_for_this_species{$genus};
			};

		if($unconstrain_broadly_labelled_data == 1)
			{
			unless($lineage =~ /family\:[A-Z][a-z]+/){$lineage = "NA"; $new_members_without_family_names_ignored++};
			};

		unless($lineage eq "NA")
			{

			# SEP2016, unconstrain these
		#	}else{
		#	die "WARNING ... cant find lineage.\n";
			
		#	print "\n$coutn of $#new_terminals, ID:$new\nlineage:$lineage\n";

			my $current_member_to_be_assigned = "NA";
			while($lineage =~ s/\:(\w+)//)
				{
				my $test_tax = $1;# print "\t\ttax in current lineage:$test_tax\n";
				if(exists($all_constraint_taxa{$test_tax}))
					{
					$current_member_to_be_assigned = $test_tax; # print "\t\t\texist as constraint taxon [re]assigning\n";
					};
				};

			$new_member_assigned_to_constraint_taxon { $new } = $current_member_to_be_assigned;
			print OUT_LIST "$new\n";
			};

	};#$new(@new_terminals)

my @testkeys = keys %new_member_assigned_to_constraint_taxon;

close OUT_LIST;

print "
constraint tax determined for each of $#testkeys new members 
new_members_without_family_names_ignored:$new_members_without_family_names_ignored

";




# finally need to read constraint tree and replace tax with new member lists
# $yet_another_newick_string



############################################################################
place_specieslevel_data_into_taxon_constraints($yet_another_newick_string);#
############################################################################


# all labels on input constrqint whic have been string replaced.:
# add_these_to_this

# all_constraint_taxa contains all seperate tax names assigned .. NOT relational constraints types (hyphenated)


$xy = 0;
if($xy == 1)
{
my @tax_on_in_constrt = keys %all_constraint_taxa;@tax_on_in_constrt = sort @tax_on_in_constrt;
print "\nlook from $#tax_on_in_constrt tax on the input constraint tree, check thses are all replaced\n";
foreach my $tax(@tax_on_in_constrt)
	{
	if(exists($add_these_to_this{$tax}))
		{
	#	print "\t$tax replaced\n";
		# taxa as node labeles will be removed at next step so dont need worry about these
		# ony worry abotu input taxa which remain which are terminals.
		}elsif($yet_another_newick_string =~ /\)[A-Za-z\-]*$tax/)
		{}elsif($yet_another_newick_string =~ /$tax[^a-zA-Z\_]/)
		{
		my $context = "NA";
		if($yet_another_newick_string =~ /(........$tax........)/)
			{$context = $1};
		print "warning. terminal taxon ($tax) in your INPUT constraint tree might be retained in output\n";
		print "\tin some cases beacuse no seqs of this tax in your fasta file.\n";
		print "\tafter removing node labels, grep both fasta file and newick, rm if neccessary\n";
		print "\tcontext in new constraint tree:$context\n";
		push @prune_tax , $tax;
		};
	};
};

# try different way::::


my @tax_on_in_constrt = keys %all_internal_node_strings_assigned_to_input_constraint;
@tax_on_in_constrt = sort @tax_on_in_constrt;
foreach my $tax(@constraints_writtenn)
	{
	if($yet_another_newick_string =~ /$tax/)
		{
		print "constraint:$tax\n";$constraints_klinging_on++;

#		if($yet_another_newick_string =~ /\)$tax/)
#			{
		#	print "\tnever mind, just node label\n";
#			}els
			if(
			# $yet_another_newick_string =~ /[^a-zA-Z\_\-\)]$tax[^a-zA-Z\_\-]/

			$yet_another_newick_string =~ /[^a-zA-Z\_\-\)]$tax[^a-zA-Z\_\-]XXXXXXXXXX/


			)
			{
			print "\tmust be terminal.... needs pruning\n";
			$tax =~ s/\-/"\\-"/e;
			push @prune_tax , $tax;
			}elsif($yet_another_newick_string =~ /$tax[\_a-zA-Z0-9]{6}/){
		#	print "\tno problem, its something else\n";
			}else{
			my $context = "NA";
			if($yet_another_newick_string =~ /(.{20}$tax.{20})/){$context = $1};
		#	print "\thuh!\n\tcontext:$context\n";
			};
		};
	};

if($constraints_klinging_on =~ /[\w\d]/)
	{
	print "\nconstraints_klinging_on:$constraints_klinging_on\n";
	};


# if( $ugly_hack == 1){ @prune_tax = ( "Dytiscoidea\\-Haliploidea" , "Cucujoidea\\-Chrysomeloidea") };


if($#prune_tax >= 0 && $process_new_consraint == 1)
	{
print "\nyou have opted to prune ... \n";

	foreach my $tax(@prune_tax)
		{
		print "terminal to prune:$tax\n";
		# context in new constraint tree:lisonia,Archaeopsyllinae,(Polype


	my $testsrting = $yet_another_newick_string;
	my $count_instances_of_tax = 0;
	while( $testsrting =~ s/[\,\(\)]$tax\,//){$count_instances_of_tax++};
	if($count_instances_of_tax >= 2){print "\nwarmning, multipple instances ($count_instances_of_tax) of tax ($tax) in newick string \$yet_another_newick_string\n"};

		if($yet_another_newick_string =~ s/\,$tax(\,)/$1/)		
			{
			print "\teasy to prune, was in a polytomie\n";
			}elsif($yet_another_newick_string =~ /\,\(($tax\,\(.+)/)
			{
			print "\tfound tax in bifurcation\n";
			# nt tree:alinae,(Chrysidini,(Poli	

			my $string = $1;
			my $count_nestedness = 1;my $string_to_replace = "(";my $regex_to_replace = "\\(";
			my $parentheses_encountered =0;

			my $x = 1;
			while($x == 1)
				{
				my $preceeding_string; my $next_parenthesis;
				if($string =~ s/^([^\(\)]*)([\(\)])//)
					{
					$preceeding_string = $1; $next_parenthesis = $2;
					}else{die "\nerror 522\n"};

				#unless($parentheses_encountered == 0){
				$regex_to_replace .= "$preceeding_string" . "\\" . "$next_parenthesis";
				$string_to_replace .= "$preceeding_string$next_parenthesis";
					#};
				$parentheses_encountered++;

				if($next_parenthesis eq "("){$count_nestedness += 1};
				if($next_parenthesis eq ")"){$count_nestedness -= 1};

			#	print "parentheses_encountered:$parentheses_encountered nestedness:$count_nestedness string:$preceeding_string next:$next_parenthesis\n";
				if($count_nestedness == 0)
					{
					$x = 0; # break loop
				#	print "break loop\n";
					};
				};


	# 	(Chrysidini,
	#		(Polistes,Abispa,Vespa,Wallacidia,
	#			(Linepithema,Leptomyrmex,Solenopsis,Camponotus,Polyrhachis,Formica,Myrmica,Pristomyrmex,Vollenhovia,
	#				(Nomada,Colletes,Osmia,Hylaeus,Melipona,Megachile,Lasioglossum,Halictus,Andrena,Seladonia,Sphecodes,Bombus
	#				)Philanthus
	#			)Vespoidea-Apoidea-Sphecoidea
	#		)


			my $replace_with = $string_to_replace;
			unless($replace_with =~ s/(\()$tax\,\((.+)\)/$1$2/){die "\nerror 555.\n"};

			my $following_string = "NA"; 
			if($yet_another_newick_string =~ /$regex_to_replace(.{40})/)
				{$following_string = $1}

		#	print "\nstring_to_replace:$string_to_replace\n";
		#	print "replace_with:$replace_with\n";
		#	print "REGEX_to_replace:$regex_to_replace\n";
		#	print "following 40 chars:$following_string\n";

			# here need to make sure non a-z0-9 characters are appropriate for use in regex,
			# for example parentheses need two backslashes before (done already)
			# maybe some hyphens . others?

			# sometimes have new tax as node label for some reason, just rm
			# us,Haliplus_fulvus,Haliplus_fasciatus)Gyrinus_sp_BOLD_AAG0707_943417329,Dineut
			$yet_another_newick_string =~ s/($regex_to_replace)[A-Z][a-z]+[A-Za-z0-9\_]+/$1/;

			if($yet_another_newick_string =~ s/$regex_to_replace/$replace_with/)#
				{
				print "\nsuccess\n";
				if($yet_another_newick_string =~ /(.{20}$tax.{20})/){my $error = $1;print "\nerror, still there!\n($error)\n"};
				}else{
				print "\nstring_to_replace:$string_to_replace\n";
				print "replace_with:$replace_with\n";
				print "within:$yet_another_newick_string\n";
				die "\nerror 560 trying to prune in newick string\n"
				};

		#	die;

			}else{
			if($yet_another_newick_string =~ /(.{20}$tax[^a-zA-Z\_].{20})/)
				{$contenxt = $1; print "\ncontext:$contenxt\n"}else{print "\ncant even find it\n"};	
			die "\nerror .. this was only a quick hack .. not written for current configuration\n"
			}

		};

	};


open (OUT67 , ">$outfile_prefix.newick_constraint_tree2") || die "\n\n";
print OUT67 "$yet_another_newick_string;\n";
close OUT67;




# do basic constraints, at the genus level .. 

my @backbone_genera = keys %terminals;@backbone_genera = sort @backbone_genera;

# used to assist pruning
# open(TEST_FILE , ">$outfile_prefix.needs_pruning") || die "\nerror 759\n";

my $j = 0;
foreach my $genus(@backbone_genera)
	{
	$j++;
	my $count_new_members_assigned_to_this_genus = 0;
	my $replace_backbone_genus_with_string = "";
	foreach my $new(@new_terminals)
		{
		if($new =~ /$genus[_]/)
			{
			$count_new_members_assigned_to_this_genus++;
			$replace_backbone_genus_with_string .= "$new,";
		#	$print "\tnew member $new\t"
			};
		};
	$replace_backbone_genus_with_string =~ s/\,$//;
	print "Genus $j of $#backbone_genera backbone genera. Setting basic constraints named $genus, " , 
		"found $count_new_members_assigned_to_this_genus\n";
	if($count_new_members_assigned_to_this_genus >= 1)
		{
		print TEST_FILE ">$genus\nACTG\n"
		}else{print "\n\tfound nothing for backbone genus $genus .. will need pruning. quitting.\n"};

	$backbone_copy_newick =~ s/$genus([^a-z])/$replace_backbone_genus_with_string$1/;

	};

close TEST_FILE;


open (OUT699 , ">$outfile_prefix.basic_constraint_tree") || die "\n\n";
print OUT699 "$backbone_copy_newick\n";
close OUT699;





print "

line 801,
FIN!
";


exit;









open(OUT57 , ">$outfile_prefix.taxon_constraints_output_LOG") || die "\nerrro\n";
foreach my $tax(@mps)
	{
	if($record_taxon_constraints_ouput_MPa{$tax}==1 && $record_taxon_constraints_ouput_MPb{$tax}==1)
		{
		print OUT57 "backbone monophyly $tax\tapplied to new matrix \n";
		}else{
		print OUT57 "backbone monophyly $tax\tnot applied to new matrix \n";
		}
	}
close OUT57;



#unless($outgroup_found == 1){die "\nerror assigning outgroup constraint.\n"}


close OUT5;

print "out of $#new_terminals new seqs, $assigned22 can be unambiguously placed on the backbone tree, 
$not_assigned22 cannot

\n";

if($no_lineage_found >= 1)
	{
	print "no lineage found for $no_lineage_found, therefor also cant assign\n";
	}else{
	print "lineage was found for all ... something\n ";
	};



%store_node_assigned_newIDs;

$newick5 = "($root_node)";
# this sub works on $newick5
#############################################
get_terminals_from_this_node30($root_node);#
#############################################

$newick5 =~ s/^\((.+)\)$/$1/;

$raxml_constraint_tree = "$treefile.rax_constraint";

open(NEWICK5 , ">$outfile_prefix.raxml_format_constraint_tree") || die "\nerror 193.\n";
print NEWICK5 "$newick5;\n";
#print "$newick5\n";
close NEWICK5;




# tnt: given up trying to get it running.

# open(TNT_FILE_HANDLE , ">$fas_file.$treefile.tnt") || die "\nerror 366\n";

%assigned_already;

# tnt format, 
#	# constraint where 2 taxa 'float'
#	force + [ a b c  (d e) ]
#tnt64 mxram 2000, rseed 1,p default-Second.nex timeout 0:20:00, echo=, taxname=, force + [17 18 19 20], constrain =,
#
#$tnt_command_string

$tnt_command_string = "tnt64 'mxram 2000, rseed 1,p tnt_input.nex timeout 0:20:00, echo=, taxname=, ";


# for fasttree constraint file (uses format similar to phylip)
$constraint_number;
%phylip_constraints;

# will be stored in hash under keys 1..$constraint_number
# ids you can get from my @all_new_IDs = keys %query_IDs;


# this builds list of tnt commands in $tnt_command_string

#############################################
get_terminals_from_this_node31($root_node);# 
#############################################



# fasttree constraint file.
# 5 2
#A        00
#B        0-
#C        10
#D        11
#E        -1



my @all_new_IDs = keys %query_IDs;@all_new_IDs = sort @all_new_IDs;

print "\nmaking constraint file for fasttree.
dimensions: taxa:" , scalar @all_new_IDs,
"
constraints:$constraint_number
\n";

# open(OUT_FT , ">$output_filename.ft") || die "\nerror 421\n";
print OUT_FT " " , scalar @all_new_IDs , " $constraint_number\n";

foreach my $id(@all_new_IDs)
	{
	print OUT_FT "$id        ";
		
#	for my $constaint_no(1..$constraint_number)
#		{
		my $state = $phylip_constraints{$id};print OUT_FT "$state";
#		}

	print OUT_FT "\n";
	};




#constrain =, outgroup 0, log tnt_logfile,  hold 2000, xmult=replication 100 drift 10 ratchet 10 hold 10, best, majority *, tchoose / ,ttags= , resample boot replications 500 from 0 , ttags, tsave *tnttrees, save*, tsave / ,export *tnt_output, log / ,quit;

$tnt_command_string .= "constrain =, log tnt_logfile,  hold 2000, xmult=replication 100 drift 10 ratchet 10 hold 10, best, majority *, tchoose / ,ttags= , resample boot replications 500 from 0 , ttags, tsave *tnttrees, save*, tsave / ,export *tnt_output, log / ,quit;'";



#print "\nrunnning tnt\n";
#system($tnt_command_string);

print TNT_FILE_HANDLE "$tnt_command_string";
close TNT_FILE_HANDLE;

print "
printed tnt
";



#$new_newick_string =~ s/\)(\w+)/ $assign_these_new_terminals_to_taxon{$1} . ")" . $1 /ge;

		if($new_newick_string =~ /(........Pierini.......)/)
			{print "dolar 1:$1\n"};


my @taxa_to_which_new_species_are_assigned = keys %assign_these_new_terminals_to_taxon;
@taxa_to_which_new_species_are_assigned = sort @taxa_to_which_new_species_are_assigned;
print "\n" , scalar @taxa_to_which_new_species_are_assigned , 
" taxon names on the backbone tree, to which new terminals will be added\n\n";

print  "new newick string:$new_newick_string;\n";


foreach my $taxon(@taxa_to_which_new_species_are_assigned)
	{
	if($taxon =~ /\w/)
	{
	my $new_terminals_for_this = $assign_these_new_terminals_to_taxon{$taxon};
	$new_terminals_for_this =~ s/\,$//;

	# hack to cut down tree for vewiing:
	$cut_down_for_viewing = 1;
	if($cut_down_for_viewing ==1)
		{
		$new_terminals_for_this =~ s/^([\_A-Za-z]+)\,.+\,([\_A-Za-z]+)$/$1 . "," . $2/e;#################
		};
	print "\ntax:($taxon)\tnew terminals:($new_terminals_for_this)\n";

	if($new_newick_string =~ /(.......[\(\)\.]$taxon[\(\)\.]......)/)
		{print "looks like:$1\n"};

	my $taxname = $taxon;
	$inculde_tax_branch_labels = 0;
unless($inculde_tax_branch_labels ==1){$taxname = ""};

	if($new_newick_string =~ s/\(\)$taxon([\,\(\)])/ "(" . $new_terminals_for_this . ")" . $taxname . $1/e)
		{
			#if($new_newick_string =~ /(........Pierini.......)/)
			#	{print "dolar 1:$1\n"};
		}elsif($new_newick_string =~ s/(\w)\)$taxon([\,\(\)])/ $1 . "," . $new_terminals_for_this . ")" . $taxname . $2/e)
			{
			}elsif($new_newick_string =~ s/\)\)$taxon([\,\(\)])/ ")," .  $new_terminals_for_this . ")" . $taxname . $1/e){#....tata))Cucujiformia)Polyp....
			}else{
#print  "\nERROR. newick string:$new_newick_string;\n";
		
		die "\nerror 266:cant find $taxon on tree\n"
			
			};

unless($inculde_tax_branch_labels ==1)
	{
#ta)))Nymphalidae,((Ant

while($new_newick_string =~ s/([\,\(\)])$taxon([\,\(\)])/$1$2/g){};

	};



	}};





open(NNS, ">$outfile_prefix.new_newick_string4");
print NNS "$new_newick_string;\n";
close NNS;

# this doesnt make sense:
# record_tree3("new_newick_string2");#
#########################





$new_newick_string3 = "($root_node3)";

###############################################
get_terminals_from_this_node3($root_node3 , 0);#
###############################################


open(NNS3, ">$outfile_prefix.new_newick_string3");
print NNS3 "$new_newick_string3;\n";
close NNS3;






die;


my @test_query_IDs = keys %query_IDs;
foreach my $id(@test_query_IDs)
	{
	#print "query:$id\n";
	unless (exists($terminals{$id}))
		{print "\n\nWARNING:member found in fasta file ($id) is absent from tree. are you sure you are using the correct pair of files?\n\n"}
	}


$terminals_belonging_to_current_node;




print "\n\n\nend of script\n";
exit;





#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################





sub record_tree2
{


my $tree1= "";
open(IN, $treefile) || die "\n\nerror 1408 cant open file $treefile\n\n";
while (my $line = <IN>)
	{
	if($line =~ /(.+\(.+)/){$tree1 = $1}
	}
close(IN);

print "\nlooking at tree:$treefile\n";

open(BCN_LOG, ">>backbone_constraints_newick_LOG") || die "";
print BCN_LOG "\nreading backbone phylogeny:$treefile\n";
close BCN_LOG;


$tree1 =~ s/ //g;

$second_backbone_copy_newick = $tree1;
$third_backbone_copy_newick = $tree1;  # for inserting barcode IDs
$fourth_backbone_copy_newick = $tree1; # for grafting source trees

# this one will have internal taxonomic node labesl:
$forth_backbone_copy_newick = $tree1;

# thus , rm supprto node labesl:
$forth_backbone_copy_newick =~ s/(\))[0-9\.]+/$1/g;

# grafting with internal labels
$backbone_copy_newick5 = $tree1;
$backbone_copy_newick5 =~ s/(\))[0-9\.]+/$1/g;




$tree_parse = 0; 	# assuming sensible tree.

if($plot_constraints == 1){$tree_parse = 1}; # user more interested in the R plot which only works with no lengths

if($tree_parse ==1)
	{
	# remove branchlengths, scientific notation, incl negative values for distance trees. example: -8.906e-05
	$tree1 =~ s/\:\-*\d+\.\d+[eE]\-\d+([\(\)\,])/$1/g; 
	# remove regular branchlengths: 0.02048
	$tree1 =~ s/\:\-*\d+\.\d+//g; 
	# remove 0 length branchlengths
	$tree1 =~ s/\:\d+//g;
	# whole number node support
	# estris)100)100,UCE_v2_Bomb    &   ris))100,UC
	while($tree1 =~ s/(\))100([\)\,])/$1$2/){};
	# 1 or 2 digits
	while($tree1 =~ s/(\))[0-9]{1,2}([\)\,])/$1$2/){};
	}else{
	print "script set to not remove all formats of branch length, this assumes sensible newick string\n";
	};


my $newick_string 	= $tree1;
my $interal_node	= 0;
$backbone_copy_newick = $tree1;






#####################################################################################

# new newick parser ....

#while ($newick_string =~ s/\(([^\(\)]+)\)(\d*)/INTERNAL_NODE_$interal_node/) # rooted bipart raxml tree, the 2 last nodes will have no boot values
#	{

while ($newick_string =~ s/\(([^\(\)]+)\)([0-9\.]*)/INTERNAL_NODE_$interal_node/) # processed sumtrees output from a load of raxml boots
	{

	my $node = $1;my $boot = $2; my $nodeID = "INTERNAL_NODE_$interal_node"; # print "nodeID:$nodeID node:$node\n";


		# seems raxml constraint needs to be unrooted
	if($newick_string =~ s/^\((INTERNAL_NODE_\d+)\:[\d\.]+(\,[A-Za-z_]+\:[\d\.]+)\)(\;)/$1$3/)
		{
		$node .= $2;
		print "\nwarning, since raxml constraint needs to be unrooted, " ,
			"tried to unroot your rooted tree. may or may not work.\n\troot node now:$node\n";
		};


	# found a pair of parentheses (with no parentheses inside).
	# this string is taken, and replaced with a single node (INTERNAL_NODE_\d)
	# node is numbered (\d in variable $interal_node), 
	# the number being incremented at the end of the loop
	# find what is at the current node (decendents), by splitting at the commas

	if($boot =~ /\d/){$nodes{$nodeID}{support} = $boot}else{$nodes{$nodeID}{support} = 100};#print "\tnodeID:$nodeID boot:$boot\n";

	my @child_nodes = split /\,/ , $node;#print "\@child_nodes:@child_nodes\n";
	$child_counts{$nodeID} = $#child_nodes;

	if($interal_node >= 6)
		{}elsif($interal_node == 5)
		{print "........\n.......\n"}else{
		print "nodeID:$nodeID\n\t$#child_nodes child_nodes:@child_nodes\n";
open(BCN_LOG, ">>backbone_constraints_newick_LOG") || die "";
print BCN_LOG "\tNEWICK nodeID:$nodeID\tcount_child_nodes$#child_nodes child_nodes:@child_nodes\n";
close BCN_LOG;
		
		};


	for $i(0 .. $#child_nodes)
		{
		#print "$child_nodes[$i]\n";
		my $bl = "NA"; 	
		if($child_nodes[$i] =~ /\:(.+)/)	# branchlength found
			{$bl = $1};
		$child_nodes[$i] =~ s/\:(.+)//;
		#print "node:$child_nodes[$i]\tbl:$bl\n ";

		# optional processing of tip labels
		unless($child_nodes[$i] =~ /INTERNAL_NODE_/)
			{
			if($remove_accessions_from_reference_species == 1)
				{$child_nodes[$i] =~ s/^([A-Z][a-z]+_[a-z]+)_.+$/$1/}
			if($process_backbone_tree_terminal_IDs == 1)
				{$child_nodes[$i] =~ s/.+_([A-Z][a-z]+_[a-z]+)$/$1/}
			if($process_backbone_tree_terminal_IDs == 2)
				{
				unless($child_nodes[$i] =~ s/^([A-Z][a-z]+)_.+$/$1/){die "\nERROR, process_backbone_tree_terminal_IDs == 2, however cant remove species from terminal:$child_nodes[$i]\n"};
				
				}
			};


		# record each decendent of the current node, in a hash nodes{ the current node }{ child number }
		$nodes{$nodeID}{$i} 			= $child_nodes[$i];
		$nodes{$child_nodes[$i]}{parent} 	= $nodeID;
		$nodes{$child_nodes[$i]}{branchlength} 	= $bl;


		# record length of branch between the current two nodes,
		# store in both directions
		my $current_child = $child_nodes[$i];
		$branchlengths{$nodeID}{$current_child} = $bl; # print "storing bl ($bl) to $nodeID / $current_child\n";
		$branchlengths{$current_child}{$nodeID} = $bl;


		# and record whether the current node is a terminal

		unless($child_nodes[$i] =~ /INTERNAL_NODE_/)
			{
			if($terminals{$child_nodes[$i]} == 1)
				{
				if($process_backbone_tree_terminal_IDs == 2){print "\nwarning, you have set script to remove species names, this assumes no conspecifics, otherwise expect crash here.\n"};
				die "\nfatal error, found duplicate tip in backbone tree:$child_nodes[$i]\n"
				};

			# terminal labels in a tree published in PNAS:
			# Papilionoidea_Nymphalidae_Nymphalinae_Polygonia_c_album_INSfrgTAJRAAPEI__20
			# Papilionoidea_Papilionidae_Papilioninae_Papilio_machaon_Papilionoidea
			if($child_nodes[$i] =~ /[A-Z][a-z]+_[A-Z][a-z]+_[A-Z][a-z]+_/){$poorly_labelled_terminals++};
			if($poorly_labelled_terminals >= 20)
				{die "\nERROR. your terminal labels (e.g. $child_nodes[$i]) cannot be parsed. should be no more than Genus_species. quitting.\n\n "};
 	
			$terminals{$child_nodes[$i]} =1;$count_the_terminal_nodes++;
			my $genus_name = $child_nodes[$i];$genus_name =~ s/^([A-Z][a-z]+)_.+/$1/;
			my $complete_lineage = $complete_lineage_for_this_species{$child_nodes[$i]};
			unless($complete_lineage =~ /\w/){$complete_lineage= $complete_lineage_for_this_species{$genus_name}};
			unless($child_nodes[$i] =~ /^[\w\_]+$/)
				{
				print "\nerror 100. no terminal label:$child_nodes[$i]\n";
				# die "";
				}
			if( $complete_lineage =~ /\w\s/)
				{
				$complete_lineages_retreived++;

				#print "child_nodes[i]:$child_nodes[$i] complete_lineage:$complete_lineage\n";
				# need to count how many terminals there are for each taxa,
				# for a given taxa, if all of them are decended from one node
				# it means they are monophlyetic, and will be constrained.

			#	print "\n$child_nodes[$i]\n";
				while($complete_lineage =~ s/^([^:]+):(\w+)//)
					{
					my $current_rank = $1;my $current_taxname = $2;
				#	print "\tcurrent_rank:($current_rank) current_taxname:($current_taxname)\n";
				#	current_rank:( family) current_taxname:(Bothrideridae)
					$count_taxa_represented_in_terminals{$current_taxname}++;
				#	$tax_is_of_this_rank
					}

				}else{
				$errors_printed++;
				if($errors_printed < 6)
				{print "\nERROR (1195). sub record_tree2. no complete lineage retrieved for genus $genus_name\n" , 
				"\t$complete_lineages_retreived have sucessfully been retrvied previously\n" , 
				"\ttip: check you are not using outdated NCBI taxonomy database; check taxon specified in taxon table ($root_taxon_name) matches those in your tree\n";
				print "check its in there with command grep \"$genus_name\" names.dmp\n";
				}elsif($errors_printed == 6){print "..... not printing further errors of this type.\n"};

				if($genus_name =~ /^Beetles$|^Flies$|^Moths$|^Wasps$/)
					{
					}else{
				#	die ""
					};
				# die "\tquitting.\n\n";
				};



			}#unless($child_nodes[$i] =~ /INTERNAL_NODE_/)

		}

	#print "node:$interal_node\n\tchild nodes:@child_nodes\n";
	$root_node = $nodeID; # this is a global and will keep the value in the final loop, so the root can be identified later.
	$interal_node++;$print_internal_nodes++;
#	if(length($newick_string)<= 500){print "newick_string:$newick_string\n"};

	}; # while ($newick_string =~ s

#####################################################################################

print "
taxonomy not retrived for $errors_printed terminals
";

unless($complete_lineages_retreived >= 1)
	{
	print "\n\n\nerror, your taxon table does not contain any of the things in your phylogeny. " , 
		"please check root tax number input to taxon_table.pl. quitting.\n\n\n";
		die "";
	};


 print "remaining of newick_string :$newick_string\n";

my @terminal_array = keys %terminals; @terminal_array = sort @terminal_array;
# store single value, used later:
$count_terminals = scalar @terminal_array;

open(BCN_LOG, ">>backbone_constraints_newick_LOG") || die "";
print BCN_LOG "finished reading phylogeny, 
  count_terminals:$count_terminals
  taxonomy not retrived for $errors_printed terminals
\n";
close BCN_LOG;


unless($interal_node >= 2){die "\nerror reading your phylogeny.\n"}


}#sub record_tree2




#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



# pseudocode:
# for each barcode taxon
#	retreive taxonomic lineage.
# 	for each taxon in lineage,
#		append barcode ID in list of its members.


sub read_fas
{

open(IN, $fas_file) || die "\nerror ... fasta file you have given ($fas_file) cannot be opened. 
make sure it is named correctly and in the working directory. quitting\n";
print "sucessfully opened fasta file $fas_file\n";


print BCN_LOG "\nopened fasta file $fas_file, reading IDs ...\n";
# close BCN_LOG;

open(OUT1111, ">lineage_assignments_to_barcodes") || die "";

my $count_diet_IDs = 0;
while (my $line = <IN>)
	{
	$line =~ s/\n//;$line =~ s/\r//;
	if($line =~ /^>(.+)/)
		{
		my $id = $1; $fas_ID_Print++; if($fas_ID_Print <= 5){print BCN_LOG "native fasta ID format:$id\n"};
#		print "\nid:$id\n";

		if($id =~ /^[a-z]+_[a-z]+$/i)
			{
			$count_binomials++;
			}elsif($id =~ /^[a-z]+_[a-z]+_[a-z]+/i)
				{
				$extended_ids++;
			#	die "\nfasta ID:$id. not writenn for non-binomials \n";
				}elsif($id =~ /^[a-z]+$/i){
				$genus_names_only++;
				};
		$index_of_fastaID{$id} = $count_diet_IDs;
		$count_diet_IDs++;
		$query_IDs{$id}= 1;
		$fasta_IDs_again{$id}= 1; # print "stored $id\n";

		# get all tax names of this species, then record +1 for each.
		my $complete_lineage = $complete_lineage_for_this_species{$id};
		unless($complete_lineage =~ /\:\w/)
			{
			my $genusname = $id;

			$genusname =~ s/AntiopalaX/Antiopala/;
			$genusname =~ s/BareaX/Barea/; # BareaX_ectadia_342779467
			$genusname =~ s/^([A-Z][a-z]+)BIOUG\d+.+/$1 . "_"/e; # ArchaeognathaBIOUG07747D02)
			$genusname =~ s/^([A-Z][a-z]+)BOLD\w+\d+.+/$1 . "_"/e; # ArchaeognathaBOLDAAG6186
			$genusname =~ s/^([A-Z][a-z]+)CCDB\d+.+/$1 . "_"/e; # BlattodeaCCDB23299D03
			$genusname =~ s/^([A-Z][a-z]+)DJM2.+/$1 . "_"/e; # BlattodeaDJM2individual3


			unless($genusname =~ s/^([A-Z][a-z]+)_.*/$1/)
				{
			#	print "\nerror trying to extrat genus name from ($id)\n";
				};
		#	print "\n no lineage found for id ($id)\ntrying genus name ($genusname)\n";

			$complete_lineage = $complete_lineage_for_this_species{$genusname};
			unless($complete_lineage =~ /\:\w/)
				{
			#	print "warning ... cant find lineage for sequence id:$id genusname:$genusname\n";
				# as of SEP2016, will try without constraints for some unlabeled BOLD data
				# i note mispellings slip though here and will be unconstrained
				# but these are mosrly the order labelled BOLD seqs
				};
			}# unless($complete_lineage =~ /\:\w/)
		if($print_ID_count <= 10)
			{
			print "fasta ID:$id\n";
			}elsif($print_ID_count == 11){print ".......\n"};
			$print_ID_count++;

	#	print "complete_lineage:$complete_lineage\n";
		if($complete_lineage =~ /\w/){$tax_lineage_found_for_fasta_entry++}else{$tax_lineage_not_found_for_fasta_entry++};

		# sometimes identical string (name) used for different ranks
		my %tax_used = ();
		print OUT1111 "$id\t$complete_lineage\n";


		while ($complete_lineage =~ s/ (\w+)\:(\w+) / /)
			{
			my $rank = $1; my $tax = $2;$count_tax_in_new_data{$tax}++; #print "tax:$tax\n";
			if($rank eq "family"){
			$family_of_fasta_ID{$id} = $tax};

			unless($tax_used{$tax} == 1)
				{
				$fasta_ids_for_each_taxon{$tax} .= "\t$id\t"; # this is prob key point at which barcode species are stored.
				};	
			$tax_used{$tax} = 1;

			};		

		my @tax_current_fasta = keys %tax_used;

		if($id =~ /^F/) {
		#	print "id:$id tax_current_fasta:$#tax_current_fasta LIN:$complete_lineage\n";
			};
		}else{
		if($line =~ /^>/){die "\nerror 646.\n"}
		}
	}
close(IN);
close OUT1111;
print BCN_LOG "\tcount_binomials:$count_binomials 
\textended_ids (Order_Genus_species):$extended_ids
\tgenus_names_only:$genus_names_only

";

if($tax_lineage_not_found_for_fasta_entry > (($tax_lineage_found_for_fasta_entry + $tax_lineage_not_found_for_fasta_entry)*0.5))
	{
	print "\nWARNING, taxonomic lineages were not found for most your fasta entries\n";
	print "  found $tax_lineage_found_for_fasta_entry out of $count_diet_IDs\n";
	print "  this can happen for instance if your fasta file is insect-wide, and your taxon table was built for a certain order\n\n";
	};







$autodetect_fasta_IDs = 0;
if($count_diet_IDs == $count_binomials)
	{$autodetect_fasta_IDs = "binomial";print "\nautodetected strict binomial used in your fasta file. good!\n"};
if($count_diet_IDs == $extended_ids)
	{$autodetect_fasta_IDs = "extended";print "\nwarning. you are not using strict binomials in your fasta file, not guarunteed to work\n"};

my @check_species_filtered = keys %query_IDs;
# print "check_species_filtered:@check_species_filtered\n";die "";

if($count_diet_IDs == scalar @check_species_filtered)
	{
	print "\npass check, all strings in fasta file are unique (species filtered)\n"
	}else{
	print "\nwarning. duplicate IDs found in your fasta file.\n"
	}

print "\n*** $count_diet_IDs fasta IDs *** in your fasta file ($fas_file). will be looking for these in the tree ...\n\n";

unless($count_diet_IDs >= 1)
	{die "\nerror ... no fasta ID's found in the file ($fas_file). quitting.\n"}





}#sub read_fas





#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



 # *** THIS SUB NOT IN USE ***


sub get_terminals_from_this_node
{
my $next = $_[0];my $sum_branchlength = $_[1];


#print "\nNEW NODE ($next)";

my $child_nodes = $child_counts{$next};
my @next1 = ();
my @non_terminal_child_nodes = ();
my @child_nodes_test = ();


my $all_terminals_from_both = "";
my @taxon_for_child_nodes = ();
my $tax_for_all_child_nodes = "";
my $all_child_nodes_are_internal =1;

my $tax_for_parent = "";
					# do this even if child node is terminal
if($next =~ /INTERNAL_NODE_\d/) 	# && $all_child_nodes_are_internal == 1)
	{
	$all_terminals_from_this_node ="";
	get_terminals_from_this_node2($next, 0);

			# no longer a sub of this name
	$tax_for_parent = get_shared_taxonomic_name($all_terminals_from_this_node);
	#print "current node ($next) is internal, shared tax name ($tax_for_parent)\n";
	};


for $i(0 .. $child_nodes)
	{
	$count_terminal_nodes = 0;
	$all_terminals_from_this_node ="";# for current child node only

	if($nodes{$next}{$i} =~ /INTERNAL_NODE_\d/)
		{
		# if current child is internal node, go and get the decending terminals
			# this sub does this for terminal nodes only:
			# $all_terminals_from_this_node .= "$test\t";
			# gets binomial names
		#print "line561. getting terminal binomials for child node $nodes{$next}{$i}\n";
		get_terminals_from_this_node2($nodes{$next}{$i}, 0);
		
		}else{
		# otherwise, just append the node name itself
		my $established_name = $nodes{$next}{$i};
		$all_terminals_from_this_node .= "$established_name\t";
		$all_child_nodes_are_internal =0
		};



		$all_terminals_from_both .= $all_terminals_from_this_node;
		$all_terminals_from_this_node =~ s/\t$//;


			# no longer a sub of this name
				##############################################################
		my $tax = 	get_shared_taxonomic_name($all_terminals_from_this_node);#
				##############################################################

		push @taxon_for_child_nodes, $tax;
		#print "child $i of node $next has ID $nodes{$next}{$i}. $count_terminal_nodestermianls got, shared tax:$tax.\n";

		# here assess if child node is same taxonomy as node.
		# if so, this node is collapsed.
		# this is done by stepping ahead to getting the two grandchild nodes, 
		# then write both these in instead of the child node.


		if($tax eq $tax_for_parent)
			{
			print "this child has same taxon ($tax,$tax_for_parent), collapsing.\n";
			my $grandchild_nodes = $child_counts{$nodes{$next}{$i}};
			for $j(0 .. $grandchild_nodes)
				{
				my $push_node = $nodes{$nodes{$next}{$i}}{$j};
				print "push grandchild to array j:$j, push_node:$push_node\n";
				push @next1 , $push_node;
				if($push_node =~ /INTERNAL_NODE_\d/)	# for makeing the newick string with termainals for new data,
					{
					push @non_terminal_child_nodes, $push_node;
					push @child_nodes_test, $push_node;
					}elsif# first need to remove existing terminals
					($push_node =~ /\w/)
						{
						push @child_nodes_test, $push_node;
						}else{die "\n\nerror 788\n"}	
				};
			}else{
			print "this child has diff taxon ($tax,$tax_for_parent)\n";
			my $push_node = $nodes{$next}{$i};
			push @next1 , $push_node;	
			if($push_node =~ /INTERNAL_NODE_\d/)
				{push @non_terminal_child_nodes, $push_node};
			push @child_nodes_test, $push_node;
			};

		if($i == $child_nodes)
			{
			$all_terminals_from_both=~ s/\t$//;
			print "line793, calling sub get_shared_taxonomic_name\n";


						# no longer a sub of this name
			$tax_for_all_child_nodes = get_shared_taxonomic_name($all_terminals_from_both);
			unless($tax =~ /\w/){die "\nno tax name retreived for this set of terminals\n"}
			};



	}




 # *** THIS SUB NOT IN USE ***




	# this is the default way to building up the newick string with new nodes.
	# it is overwritten later in situations where the node needs to be collapsed.


# this is collapsed where appropriate, but includes terminal nodes. 
# printed to an additional newick string so user can see better how tree has been collapsed 
my $join_child_nodes_test = join ',', @child_nodes_test;

# removed terminals:
my $join_child_nodes = join ',',@non_terminal_child_nodes;
my $swap_string ="";

#if($#non_terminal_child_nodes <= 0)
#	{
#	print "no decendents. dead end.\n";
##$swap_string = "$tax_for_parent";
#$swap_string = "($join_child_nodes)$tax_for_parent";
#	}else{
$swap_string = "($join_child_nodes)$tax_for_parent";			# removed terminals
my $swap_string_test = "($join_child_nodes_test)$tax_for_parent";	# no removed terminals
#	}


#print "tax_for_parent:$tax_for_parent\n";
#die;
$all_internal_taxa_on_constraint_tree{$tax_for_parent}=1;


if($next =~ /INTERNAL_NODE_\d/ && $all_child_nodes_are_internal == 1)
	{
#	$all_terminals_from_this_node ="";
#	get_terminals_from_this_node2($next, 0);
#	$tax_for_parent = get_shared_taxonomic_name($all_terminals_from_this_node);


	my $all_branches_from_this_node_are_the_same_taxon = 1;
	foreach my $test_taxon(@taxon_for_child_nodes)
		{
		if($tax_for_parent eq $test_taxon)
			{}else{$all_branches_from_this_node_are_the_same_taxon=0}
		};


	# collapse node method 1
	# turns not really what i want to do here.
	# but i keep the code in case it is what i want to do in a future situation
	if($all_branches_from_this_node_are_the_same_taxon ==1)
		{
		print "node $next will be COLLAPSED. tax_for_parent:$tax_for_parent taxon_for_child_nodes:@taxon_for_child_nodes\n";
	#	$swap_string = "$join_child_nodes";
		# collapse node. get identity of parent of current node (which is itself parent)
		# then add the child nodes to this  
		# no. instead you need to get the grandchild node of the matching taxon branch, and pull it up to be a child node.

		#$child_counts{$nodeID} = $#child_nodes;
		#$nodes{$child_nodes[$i]}{parent} 	= $nodeID;

		}


	};

my $join_the_child_nodes;
$collapse_nodes = 0;
if($collapse_nodes == 1)
	{
	# build the newick string, replace current internal node with child nodes
	print "node:$next becomes ($swap_string)\n";
	$new_newick_string =~ s/$next(\W)/$swap_string$1/;		# has terminals removed
	$newick_print_string=~ s/$next(\W)/$swap_string_test$1/;
	}else{
	@next1 = ();
	for $i2(0 .. $child_nodes)
		{push @next1, $nodes{$next}{$i2}};
	# default:
	$join_the_child_nodes = join ',', @next1;
	$swap_string = "($join_the_child_nodes)$tax_for_parent";
	$new_newick_string =~ s/$next(\W)/$swap_string$1/;
	$newick_print_string=~ s/$next(\W)/$swap_string$1/;

	};

# collapsing is acheived by not looping to child nodes, but skipping strait to grandchild nodes.
# similarly, grandchild nodes are printed to the growing newick string instead of child nodes.


for my $index(0 .. $#next1)
	{
############	my $test = @next1[$index];#print "test:$test\n";####################################### surely this is wrong?
	my $test = $next1[$index];

	if($test =~ /^INTERNAL_NODE_/)
		{
		#####################################
		get_terminals_from_this_node($test , $sum_branchlength+$branchlength_to_child );#	# recurse
		#####################################
		}
	}

return();

 # *** THIS SUB NOT IN USE ***



}






#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################






sub get_shared_taxa
{
my $tax_list = shift;
my @tax_array = split /\t/ , $tax_list;
my $shared_tax_substring = "";

print GRAFT_LOG "\t\tcomplete list of taxa fed into sub get_shared_taxa:@tax_array\n";

# gets lineage (series of taxa) shared by all members. also option to get only lowest shared taxa


####################################################################################################################
# 2022-08-06: previous implementation took the first terminal of list, inferred lineage to be compared to remainder.
#		this assumes near complete taxonomic inference for terminals 
#		(since if lineage not retreived for first in list, process would break),
#		which is not the case for parsing source trees which are regulaly non-molecular.
#		thus need to instead look for first instance of lineage across terminal list.
#
my $test_taxonomy = "";
my $terminal_list_index =0;
for my $index(0 .. $#tax_array)
	{
	$terminal_list_index = $index; # print "\nsub get_shared_taxa, $index OF $#tax_array\n";
	my $test_species = $tax_array[$index];

	# hash complete_lineage_for_this_species, has not just species, 
	# but has lineage for every taxon in ncbi, only modification is space->underscore
	$test_taxonomy = $complete_lineage_for_this_species{$test_species};
	#print "line928. test_species:$test_species test_taxonomy:$test_taxonomy\n";

	unless($test_taxonomy =~ /\w+/)
		{
		my $genusname = $test_species;$genusname =~ s/^([A-Z][a-z]+)_.+/$1/;
		$test_taxonomy = $complete_lineage_for_this_species{$genusname};
		# print "cant find taxonomy for species ($test_species) trying genus ($genusname)\n";
		if($test_taxonomy =~ /\w+/)
			{
			# print "\tfound using genus name ($genusname), continuing ...\n";
			}else{
			if($test_species =~ /optera_[A-Z]/)
				{print "\nlooks like user format was not processed to an established tax name\n"};
		#	 print "\nerror 1701, no taxonomy found for species ($test_species) nor genus ($genusname). test_species:($test_species) tax_list:($tax_list) \nquitting ...\n";
			$skipnode=1;
			}
		}
	if($test_taxonomy =~ /\w+/)
		{
	#	print "found first lineage at terminal index $index\n";
		last
		}else{
	#	print "keep looking\n"
		};
	};
####################################################################################################################

print GRAFT_LOG "\t\tfirst lineage retreival in list, will compare taxa of this to remaining terminals. Terminal index $terminal_list_index, Linage:$test_taxonomy\n";

my $most_inclusive_name = "NA";
my $most_inclusive_lineage = "";


if($terminal_list_index == $#tax_array)
	{
	# only a single lineage will be found for this list of terminals,
	# which cannot really be called a shared taxon, so return nothing
	print GRAFT_LOG "\t\tonly a single lineage will be found for this list of terminals, which cannot really be called a shared taxon, so return nothing\n";
	return($most_inclusive_name);
	};

my $skipnode=0;

if($skipnode==0)
{

while($test_taxonomy =~ s/^([^:]+):(\w+)//)
	{
	my $current_rank = $1;my $current_taxname = $2;
	# over counts, should not count taxa at parent node:
	#$number_of_nodes_to_which_this_taxon_has_been_assigned{$current_taxname}++;
	#print "\tcurrent_rank:$current_rank current_taxname:$current_taxname\n";

	if($#tax_array >= 1)
		{
		# this is the regular situation, there is more than one terminal you have got a lineage for,
		# and you are finding the shared taxa for these

		my $all_members_have_this_name =1;

#		foreach my $tax(@tax_array)
		for my $terminal_index($terminal_list_index .. $#tax_array)
			{
			my $tax = $tax_array[$terminal_index];
			my $test_taxonomy2 = $complete_lineage_for_this_species{$tax};
			my $genus_name = $tax;$genus_name =~ s/^([A-Z][a-z]+)_.+/$1/;
			unless($test_taxonomy2 =~ /\w/){$test_taxonomy2 = $complete_lineage_for_this_species{$genus_name}};
			#print "\tspecies:$tax complete taxonomy:$test_taxonomy2\n";

			# for current taxonomic name in lineage list, see if another terminal also has this
			if($test_taxonomy2 =~ /\w/)
			{
			if($test_taxonomy2 =~ /\:$current_taxname\s/)
				{
				print ""
				}else{
				$all_members_have_this_name =0;#print "0";
			#	print GRAFT_LOG "\t\t\tterminal $tax is not $current_taxname, its lineage:$test_taxonomy2\n";
				}
			};
			}

		#print "\nall_members_have_this_name:$all_members_have_this_name\n";

		if($all_members_have_this_name == 1)
			{
			$most_inclusive_name = $current_taxname;
			$most_inclusive_lineage .= "$current_taxname\t";
			};

		print GRAFT_LOG "\t\tcurrent_rank:$current_rank current_taxname:$current_taxname all_members_have_this_name:$all_members_have_this_name\n";

		}else{

		# alternativly, you only have one taxa, so just parse the whole lineage for this		

		$most_inclusive_name = $current_taxname;
		$most_inclusive_lineage .= "$current_taxname\t";

		}

	#print "\tmost_inclusive_lineage:$most_inclusive_lineage\n";
	}


};

# print name only:
return($most_inclusive_name);
# or lineage:
#print "\treturning:($most_inclusive_lineage)\n";
# return($most_inclusive_lineage);

}





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

print "
\n\n
     ***** relational_constraints.pl *****
		  
\n";




if($arguments =~ /-treefile\s+(\S+)/)
	{
	$treefile = $1;
	}else{
	print "\nerror reading command arguments  (-treefile)\n\n";die""
	}

if($arguments =~ /-outfile_prefix\s+(\S+)/)
	{
	$outfile_prefix = $1;
	}else{
	print "\nerror reading command arguments  (-outfile_prefix)\n\n";die""
	}

if($arguments =~ /-taxon_table\s+(\S+)/)
	{
	$taxon_table = $1;
	}else{
	print "\nerror reading command arguments  (-taxon_table)\n\n";die""
	}

if($arguments =~ /-graft_source_trees_START\s+(.+)\s+graft_source_trees_END/)
	{
	$graft_source_trees_list = $1;
	}elsif($arguments =~ /graft_source_trees/)
	{
	die "\ncommand error, graft_source_trees option indicated but cant read source tree list\n";
	};

if($arguments =~ /-support\s+(\S+)/)
	{
	$support_cutoff = $1;
	}else{
	$support_cutoff = 50;
	print "user did not give cutoff for acceptable boot support (not relevent if you are not using bootstrapped trees)\nusing default ($support_cutoff).\n"
	}

if($arguments =~ /-seqfile\s+(\S+)/)
	{
	$fas_file = $1;
	}else{
	}

if($arguments =~ /-retain_terminal\s+(\S+)/)
	{
	$retain_terminal = $1;
	}else{
	}

if($arguments =~ /-plot_constraints/)
	{
	$plot_constraints = 1;
	};

if($arguments =~ /-retain_backbone_node_labels/) # applies only to grafting algorithm
	{
	$retain_backbone_node_labels = 1;
	};

if($arguments =~ /-references\s+(\S+)/)
	{
	$reference_file = $1;
	}else{
	}

if($arguments =~ /-backbone_terminal_format\s+(\d)/)
	{
	$process_backbone_tree_terminal_IDs = $1;

# $process_backbone_tree_terminal_IDs = 2;


	}else{
	die "\n\ncommand error, please give format you have used for terminal labels of backbone tree\n" , 
	"  -backbone_terminal_format [0,1,2]\n" , 
	"  if == 0, then just genus names\n", 
	"  if == 1 (unlikely), then process ids in backbone tree, e.g. Coleoptera_Apatides_fortis to Apatides_fortis\n", 
	"  if == 2, then Genus_speciesID (which is converted to genus)\n\n";

	}


if($arguments =~ /-outgroup\s+(\S+)/)
	{
	$outgroup = $1;
	unless($outgroup =~ /[A-Z][a-z]+/){die "\noutgroup switch used but cannot parse name \n"}
	}else{
	print "\nyou did not specify the outgroup.\n";
	}


if($arguments =~ /-constrain_ranks\s+(\S[^\-]+)/)
	{
	$constrain_ranks1 = $1;
	@constrain_ranks = split /\s+/ , $constrain_ranks1;
#	print "\nuser specified ranks at which to make constraints:\n@constrain_ranks\n";
	}else{
	@constrain_ranks = ("order","suborder","infraorder","superfamily","family","subfamily","genus");
#	print "\nyou did not specify ranks at which to make constraints, applying default:\n@constrain_ranks\n";
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





sub traverse_nodes
{
my $current_node = $_[0];my $current_node_taxstring = $_[1];# $current_node_taxstring to be deprecated



 # for the current node, read the list of child nodes
my $child_nodes = $ncbi_nodes{$current_node}{child_nodes};$child_nodes =~ s/^\t//;
my @child_nodes_array = split(/\t/, $child_nodes);
$nodes_traversed++;


 # seems mostly discontinued:
if($current_node == $starting_node)
	{
	my $rank = $ncbi_nodes{$starting_node}{rank};$rank =~ s/\s/_/g;$how_many{$rank}++;
#	my $current_node_taxstring = substr $current_node_taxstring, 0,1;
#	my $originalname = $ncbi_nodes{$starting_node}{name};$originalname =~ s/[\s\t]/_/g;
#	my $child_complete_lineage = $ncbi_nodes{$starting_node}{complete_lineage} . "$rank:$originalname ";
	};



foreach my $child(@child_nodes_array)
	{
	unless($current_node == $child) # one of the child nodes of the root node (1), is also 1 (NCBI system)
	{
	my $rank = $ncbi_nodes{$child}{rank};$rank =~ s/\s/_/g;$how_many{$rank}++;
	my $name_string = $ncbi_nodes{$child}{name};$name_string =~ s/\s+/_/g; ####################### sep2013 .. fixed mar2015
	my $rankcode = $ncbi_nodes{$child}{rank_code};$rank_codes{$name_string} = $rankcode;


	unless($ncbi_nodes{$current_node}{complete_lineage} =~ /\w/)#22oct2015:otherwise where shared taxon is same as basal node, wont be assigned 
		{$ncbi_nodes{$current_node}{complete_lineage}="root:$starting_name_ff "};

	my $child_complete_lineage = $ncbi_nodes{$current_node}{complete_lineage} . "$rank:$name_string ";#prev assigned $name


	# NOV 2018: user decides which ranks are constrained	
	my $constrain_this_node = 0; # print "\n\n";
	foreach my $user_constrain_rank ( @constrain_ranks )
		{
		if($rank eq $user_constrain_rank)
			{$constrain_this_node = 1};
	#	print "rank:$rank user:$user_constrain_rank constrain:$constrain_this_node\n";
		};
#	if($constrain_this_node == 1)
#		{
		$ncbi_nodes{$child}{complete_lineage} = $child_complete_lineage;
#		};

# forgot about this, explains some imperfect species level trees
#	unless($rank eq "order") # 31 aug 2016: way to stop seqs only given order ident from being constrained.
#		{
#		$ncbi_nodes{$child}{complete_lineage} = $child_complete_lineage;
#		};


	#print "storing complete lineage for taxon:$name_string\n";	
	$complete_lineage_for_this_species{$name_string}=$child_complete_lineage;
#	print "child_complete_lineage:$child_complete_lineage\n";

	$ncbi_tax_number_for_this_species{$name_string}=$child;

#	print "complete lineage:$ncbi_nodes{$child}{complete_lineage}\n";

	my $originalname = $ncbi_nodes{$child}{name};$originalname =~ s/[\s\t]/_/g;
	$ncbi_taxnumber_for_taxname{$originalname} = $child;
#	print "ncbi_number:$child tax name:$originalname\n";

#	if ( $rank eq "species" ){};

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




sub get_terminals_from_this_node2
{
my $next = shift;

my $child_nodes = $child_counts{$next};
my @next1 = ();
for $i(0 .. $child_nodes)
	{
	push @next1 , $nodes{$next}{$i};	
	#print "line1423node:$next i:$i child:$nodes{$next}{$i}\n";
	}

# sourcetree_nodes

for my $index(0 .. $#next1)
	{
	my $test = @next1[$index];#print "test 1430:$test\n";

	if($test =~ /^INTERNAL_NODE_/)
		{

		#####################################
		get_terminals_from_this_node2($test );#	# recurse
		#####################################

		}else{
		$count_terminal_nodes++;

		if($test eq "")
			{
			print "\nnext:($next)\n";
			die "\nerror 1443. looks like sub get_terminals_from_this_node2 was called with a terminal node.\n"
			}


		### here set the format used in the previous (mtgenome) tree.
		# currently im using Order_Genus in the tree, and parsing the genus name.
#		unless($test =~ s/.+_//)
#			{die "strange format:($test)\n"}

		# 8mar, Coleoptera_Tetraphalerus_bruchi and parsing genus
	#	unless($test =~ s/.+_(.+_.+)/$1/)
	#		{die "error 1229. expecting Order_Genus_species, got strange format:($test)\n"}
		# now process ids in intial reading of the tree, so they just have a regular taxon name

		$all_terminals_from_this_node .= "$test\t";

		#print "line1459. test:$test\t";
		}
	}

return();

}























sub record_tree3  # sub not used
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


sub get_terminals_from_this_node3
{
my $next = $_[0];my $sum_branchlength = $_[1];

my $child_nodes = $child_counts3{$next};
my @next1 = ();

for $i(0 .. $child_nodes)
	{
	$count_terminal_nodes = 0;
	my $push_node = $nodes3{$next}{$i};
	push @next1 , $push_node;	
	}

# default:
my $join_child_nodes = join ',', @next1;
$swap_string = "($join_child_nodes)$tax_for_parent";
$new_newick_string3 =~ s/$next(\W)/$swap_string$1/;

for my $index(0 .. $#next1)
	{
	my $test = @next1[$index];#print "test:$test\n";

	if($test =~ /^INTERNAL_NODE_/)
		{
		#####################################
		get_terminals_from_this_node3($test , $sum_branchlength+$branchlength_to_child );#	# recurse
		#####################################
		}
	}

return();

}








#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




sub traverse_backbone_tree_and_infer_taxonomic_node_labels 	# 	UNUSED HERE (but is used in taxon_constraints)
{

my $next = $_[0];my $sum_branchlength = $_[1];


print "\nNEW NODE ($next)\n";

my @non_terminal_child_nodes = ();
my @child_nodes_test = ();
my $all_terminals_from_both = "";
my @taxon_for_child_nodes = ();

my $lowest_taxname_for_node = "";
my $lineage_for_node = "";
$all_terminals_from_this_node ="";
$taxa_for_terminals_of_node = "";

# need a list of taxa derived from currnet node. ($all_terminals_from_this_node .= "$test\t";)
# get tax names for these.
# then go through all remaining members of the tree, (%terminals, not incl $all_terminals_from_this_node)
# if none of these have these tax names (get_shared_taxa), they are defined as a constraint.
# then all new members from these, are constrained.			


		# do this even if child node is terminal
if($next =~ /INTERNAL_NODE_\d/) 	# && $all_child_nodes_are_internal == 1)
	{

	# traversing backbone tree, sub run for current node.
	# get members (terminals) of this node

	##########################################
	get_terminals_from_this_node2($next, 0);#
	##########################################


	# following gets lineage (series of taxa) shared by all members of clade. 
	# also option to get only lowest shared taxa

				####################################################
	$lineage_for_node = 	get_shared_taxa($all_terminals_from_this_node);#
				####################################################
	unless($lineage_for_node =~ /\w/){die "error 1944 could not parse lineage for node ($next)\n"};
	#print "get_all_taxa 1\n";

					###################################################
	$taxa_for_terminals_of_node = get_all_taxa($all_terminals_from_this_node);#
					###################################################
	unless($taxa_for_terminals_of_node =~ /\w/){die "error 1950 could not parse taxa for terminals of node ($next)\n"};


	#########################################################
	my $monophyly_name = test_for_monophylies($all_terminals_from_this_node);#
	#########################################################

	# if monophyly found, record the taxonomic name on the backbone node.
	unless($monophyly_name eq "NA")
		{$node_monophylys{$next} = $monophyly_name};

	}else{#if($next =~ /INTERNAL_NODE_\d/)

	# for tips, you still need to know the lineage. but it wont be a shared lineage.
	# in this case you dont have a list of terminals, theres just one, the curent node

				#########################
	$lineage_for_node = 	get_shared_taxa($next);#
				#########################
	unless($lineage_for_node =~ /\w/){die "error 1963 could not parse lineage for node ($next)\n"};

	#print "get_all_taxa 2\n";
					###################################################
	$taxa_for_terminals_of_node = get_all_taxa($next);# bug fix 22 oct 2015
					###################################################
	unless($taxa_for_terminals_of_node =~ /\w/){die "error 1969 .could not parse taxa for terminals of node ($next)\n"};


	};

# parse all names, put into this object
my $lineage_string = $lineage_for_node;
while($lineage_string =~ s/([^\t]+)\t$//)
	{
	my $taxon_for_node = $1;$all_taxa_on_backbone_tree{$taxon_for_node}=1;
	#print "taxname1:$taxon_for_node lineage_string:$lineage_string\n";
	};


# for current node, parse lowest shared taxa from the lineage:
$lowest_taxname_for_node = $lineage_for_node;
unless($lowest_taxname_for_node =~ s/.+\t([^\t]+)\t$/$1/)
	{
	unless($lowest_taxname_for_node =~ s/^$starting_name_ff	/$starting_name_ff/ )# basal TOL will have shared taxon being cellular orgs, same as root specified by user
		{
		print "next:$next\nlineage_for_node:$lineage_for_node\ntaxa_for_terminals_of_node:$taxa_for_terminals_of_node\nerror 1982\n";
		die "\neror 1550, lowest_taxname_for_node:($lowest_taxname_for_node)\n"
		};
	};

# check for error:
unless($lowest_taxname_for_node =~ /^[A-Za-z\_]+$/)
	{die "\nerror 1638. lowest shared taxa for current node ($next) has unexpected format:($lowest_taxname_for_node).\n"}
# and store the LOWEST name in this hash (no intermediates)
$node_label{$next} = $lowest_taxname_for_node;

print "lowest_taxname_for_node:$lowest_taxname_for_node\n";

$all_taxnames_for_terminals_decended_from_this_node{$next} = $taxa_for_terminals_of_node;

# 
# next compare the complete shared lineage of the current node, with that of the parent.
# in the example below, the node in which all decendents are Lucanidae,
# in addition to being the position at which new Lucanids are added
# can also accomodate the above ranks Scarabaeoidea, Polyphaga,
# but not Coleoptera, which is best assigned to the parent node.
# therefor find (if present) all intermediate ranks.
#        
#                    -------
#                    |
#           -----Lucanidae
#           |        |
#           |        ------
# ---Coleoptera      
#           |        -----
#           |        |
#           -----Adephaga
#                    |
#                    -----


my $parent_node = $nodes{$next}{parent};
my $parent_taxon = "";
my $all_potential_taxa ="";

if($node_label{$parent_node} =~ /\w/)
	{
	$parent_taxon = $node_label{$parent_node};
	#print "node($next), tax($lowest_taxname_for_node) parent($parent_node) Tax($parent_taxon)\n";


	if($lowest_taxname_for_node eq $parent_taxon)	# compare the lowest tax name of current node
		{						# to lowest tax name of parent node

		# Node and parent both have the same assigned taxon, nothing to do.
		# Although whilst here, need to count the number of nodes each taxon has been assigned to 
		$number_of_nodes_to_which_this_taxon_has_been_assigned{$lowest_taxname_for_node}++;

		}else{
		while($lineage_for_node =~ s/([^\t]+)\t$//)
			{
			# go through whole lineage of current node 
			# (including lowest, which is not trimmed from this particular variable), 
			# and stop when reaching the taxon assigned to parent
			# this direction:Corcyra_cephalonica->Corcyra->Galleriinae->Pyralidae

			my $intermediate_taxon_for_node = $1;
			unless($intermediate_taxon_for_node =~ /^[A-Za-z\_]+$/){die "\nerror 1691. unexpected taxon name:$intermediate_taxon_for_node\n"}
			#print "\ttaxname1:$intermediate_taxon_for_node. comparing to that assigned to parent ($parent_taxon)\n";

			# the first element is the lowest taxon of current node
			# if that is the same as that of the parent node, 
			# there will be nothing put in $all_potential_taxa
			# if there is only one item in $all_potential_taxa, it will be same as lowest taxon
			if($intermediate_taxon_for_node eq $parent_taxon)
				{
				#print "same, break.\n";
				$lineage_for_node = "";# delete lineage string to break loop.
				$all_potential_taxa =~ s/[\s\.]+$//;# if names have been assigned, there will be terminal char
				}else{

				# taxon assigned to parent node has not been reached yet, so append taxon name 
				$all_potential_taxa .= "$intermediate_taxon_for_node.";#print "appending\n";

				# And whilst here, keep counting the number of nodes each taxon has been assigned to 
				$number_of_nodes_to_which_this_taxon_has_been_assigned{$intermediate_taxon_for_node}++;
				};
			};
		}

	}else{#if($node_label{$parent_node} =~ /\w/)
	
	unless($next eq $root_node)# only the root node wont have a parent, so no taxon will be found for that one.
		{
		die "\nerror 1703. no taxon assigned to parent ($parent_node) of node ($next)\n";
		};
	}


my $label_string_to_assign = "";
if($all_potential_taxa =~ /\w[\.\s]\w/)# taxa found which are intermediate between those assigned to node and parent
	{
	$label_string_to_assign = $all_potential_taxa	# $all_potential_taxa has fullstop seperated taxon names
	}else{							# so they can be printed by tree viewers.

	$label_string_to_assign = $lowest_taxname_for_node # otherwise, just use lowest shared taxon name
	};
unless($label_string_to_assign =~ /\w\w+/){die "\nwhy no taxon assigned to current node ($next). \n"}


my $child_nodes = $child_counts{$next};
my @next1 = ();

# record the taxonomic name(s) given to this node
$nodes{$next}{node_taxonomic_label}=$label_string_to_assign;

#print "\n\nNEW BACKBONE NODE ($next) assigned name(s):$label_string_to_assign\n";







for $i(0 .. $child_nodes)
	{
	# specify new array of the child nodes as normal for recursing
	my $push_node = $nodes{$next}{$i};
	push @next1 , $push_node;	
	}

#print "node:($next) child nodes:(@next1) label_string_to_assign:($label_string_to_assign) lineage_for_node:($lineage_for_node)\n";

#if($next eq "Apis_florea"){die""};

# get to the tip,
#NEW NODE (Apis_florea)
#node:(Apis_florea) child nodes:() label_string_to_assign:()

# not used any more ?:
# build newick string. this string is an intermediate, 
# it contains empty parentheses ready to be filled with the species-dense data,
#thus not readable by tree viewers

my $swap_string 	= "";
$join_the_child_nodes 	= join ',', @next1;
$swap_string 		= "($join_the_child_nodes)$label_string_to_assign";
$new_newick_string 	=~ s/$next(\W)/$swap_string$1/;

# alternative newick string that is readable at this stage:

if(exists($child_counts{$next}))
	{
	$new_newick_string2 	=~ s/$next(\W)/$swap_string$1/;
	}else{
	# instead of replacing node (which is now a species name) with child nodes (which there arnt any)
	# replace node with lineage
	$new_newick_string2 	=~ s/$next(\W)/$label_string_to_assign$1/;
	}


$all_internal_taxa_on_constraint_tree{$tax_for_parent}	= 1;

# the approach in this sub, is called all the way up to tip nodes, so tax lineages etc can be easier assessed.
# which means it needs to determine the end has been reached below (no child nodes for this tip)
# its a bit different from the usual approach, in which the sub is not called for tip nodes.

if(exists($child_counts{$next}))
	{

	for my $index(0 .. $#next1)
		{
		my $test = $next1[$index];

		#####################################
		traverse_backbone_tree_and_infer_taxonomic_node_labels($test , $sum_branchlength+$branchlength_to_child );#	# recurse
		#####################################

		}

	};

return();

}






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



sub get_terminals_from_this_node30	#	CURRENTLY UNUSED
{
my $next = $_[0];
my $child_nodes = $child_counts{$next};

#print "\nNEW NODE:($next)\n";

my @next1 = ();


my $add_to_current_node = "";
my $check1=0;
my $backbone_child_nodes_absent_from_new_IDs=0;
my $count_new_IDs_added_to_current_node = 0;

for $i(0 .. $child_nodes)
	{


	# get child nodes as normal
	my $push_node = $nodes{$next}{$i};

	unless(exists($query_IDs{$push_node}))
		{
		unless($push_node =~ /INTERNAL_NODE_/)
			{
			$backbone_child_nodes_absent_from_new_IDs++;
			print "warning 1875. node $push_node is not internal, yet is not found in list of new IDs\n";
			};
		}

	push @next1 , $push_node;	
	#print "child node:$i node ID:$push_node\n";

	my $node_taxonomic_label = $nodes{$push_node}{node_taxonomic_label};

	my @split_node_tax_labels = split /\./ , $node_taxonomic_label;
	foreach my $candidate_tax(@split_node_tax_labels)
		{
		#print "\tchild node is assigned candidate_tax:$candidate_tax\n";# Endopterygota->Acyrthosiphon_pisum->Acyrthosiphon

		unless($candidate_tax =~ /^[A-Za-z\_]+$/){die "\nerror 1864. unexpected taxon name:$candidate_tax\n"}

		if(exists($assign_these_new_terminals_to_taxon{$candidate_tax}))
			{
			my $string_length_new = length($assign_these_new_terminals_to_taxon{$candidate_tax});

			#print "child $i of node $next has tax label:$candidate_tax\n";
			#print "\t\tthis taxon has new IDs to be added, length:$string_length_new\n";

			# these are new IDs, comma seperated. will have comma terminating string.
			my $new_members = $assign_these_new_terminals_to_taxon{$candidate_tax};# . ",";
			$new_members =~ s/\,$//;

			# each new member come across, count in hash.
			# later check duplicates are not met.
			my @split_new_members= split /\,/ , $new_members;
			foreach my $newID(@split_new_members)
				{
				$store_node_assigned_newIDs{$newID}++;$count_new_IDs_added_to_current_node++;
				unless($newID eq $push_node){$add_to_current_node .= $newID . ",";}
				};

			# append, this will be a new child of the current node
			#$add_to_current_node .= $new_members . ",";
			my @check_split= split /\,/ , $add_to_current_node;

			#if($add_to_current_node =~ /pisumAcyr/){die "\n1887.\n"}

			#$check1++;

			}#if(exists($assign_these_new_terminals_to_taxon{$candidate_tax}))

		};#foreach my $candidate_tax(@split_node_tax_labels)

	};

#print "
#count_new_IDs_added_to_current_node:$count_new_IDs_added_to_current_node
#backbone_child_nodes_absent_from_new_IDs:$backbone_child_nodes_absent_from_new_IDs
#";
# check any havnt been assigned elsewhere (shouldnt have)
my @all_new_IDs_assigned_at_this_point = keys %store_node_assigned_newIDs; 
foreach my $newID(@all_new_IDs_assigned_at_this_point)
	{
	if($store_node_assigned_newIDs{$newID}>1)
		{
		print "\nerror 1880. new ID:($newID), to be assigned to node:($next),\n";
		print "which has tax:(...), was assigned previously. quit.\n";die};
	#unless($newID eq ""){};# terminal comma needs keeping
				# so last object will be empty
	};




my $copy = $add_to_current_node;#$copy .= join ',', @next1;
#$copy =~ s/\,$//;
my @check_split= split /\,/ , $copy;
my %check_hash=();
foreach my $check4(@check_split)
	{
#unless($check4 =~ /^[A-Z][a-z]+_[a-z]+$/){die "\nerror 1828:$check4\n"}
	if(exists($query_IDs{$check4}))
		{
	$check_hash{$check4}++;#if($check_hash{$check4}>= 2){die "\nerror 1822, $check4 assigned already\n"}
	#$check_hash4{$check4}++;#if($check_hash4{$check4}>= 2){die "\nerror 1824, $check4 assigned already\n"}
		}
	}

my @check_split2= keys %check_hash , @next1;
my %check_hash2=();
foreach my $check4(@check_split2)
	{

	$check_hash2{$check4}++;#if($check_hash4{$check4}>= 2){die "\nerror 1824, $check4 assigned already\n"}

	}


my @new_nextnode = keys %check_hash2;

#if($check1>= 3){die "\nerror 1816\n"}


#print "node:$next child nodes:@next1 node_taxonomic_label:$node_taxonomic_label\n";
if($add_to_current_node =~ /./)
	{
	#print "new members will be added to this node, ($#check_split)\n"
	};
#if($next1[1] =~ /Trachypachus_holmbergi/i){die}

#	


$howmany+= $backbone_child_nodes_absent_from_new_IDs;

if(
$backbone_child_nodes_absent_from_new_IDs >= 1 && 
$count_new_IDs_added_to_current_node < 2
){
	# this adds to 2, so ok
unless($count_new_IDs_added_to_current_node==1 && $backbone_child_nodes_absent_from_new_IDs == 1)
	{
unless($count_new_IDs_added_to_current_node == 0 && $backbone_child_nodes_absent_from_new_IDs == 1)
{
#$howmany++;
#die "\nyou have a problem. current node ... \n";
}	}
}

# default:
my $join_child_nodes = $add_to_current_node . join ',', @next1;

$swap_string = "($join_child_nodes)";
unless($newick5 =~ s/$next(\W)/$swap_string$1/)
	{print "\nerror 1962, cant build string\n"};

for my $index(0 .. $#next1)
	{
	my $test = $next1[$index];#print "test:$test\n";

	if($test =~ /^INTERNAL_NODE_/)
		{
		#####################################
		get_terminals_from_this_node30($test  );#	# recurse
		#####################################
		}
	}

return();

}








#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



sub get_terminals_from_this_node31 	#  CURRENTLY UNUSED
{
my $next 	= $_[0];
my $child_nodes = $child_counts{$next};


my @next1 = ();
my $add_to_current_node = "";
my $check1=0;
my $backbone_child_nodes_absent_from_new_IDs=0;
my $count_new_IDs_added_to_current_node = 0;

my $node_taxonomic_label = $nodes{$next}{node_taxonomic_label};

# keep this print: 
print "\n*** NEW NODE of backbone *** 
Making TNT / FT constraint no:$constraint_number; node:($next); tax:($node_taxonomic_label)\n";


 # sometimes they have multiple names a la: Diabrotica.Diabroticites.Diabroticina.Luperini.Galerucinae
my @split_node_tax_labels = split /\./ , $node_taxonomic_label;

my $new_IDs_going_to_this_node = "";
my $count_tax_assigned_to_node =0;
my $total_tax_assigned_to_node  =scalar @split_node_tax_labels;
foreach my $candidate_tax(@split_node_tax_labels)
	{
	$count_tax_assigned_to_node++;
	print "\ttax:$count_tax_assigned_to_node of total:$total_tax_assigned_to_node assigned. taxon name:$candidate_tax\n";# Endopterygota->Acyrthosiphon_pisum->Acyrthosiphon

	unless($candidate_tax =~ /^[A-Za-z\_]+$/){die "\nerror 1864. unexpected taxon name:$candidate_tax\n"}

	if(exists($assign_these_new_terminals_to_taxon_TNT{$candidate_tax}))
		{
		print "\tIS new members to assign to this name\n";

		if(exists($assigned_already{$candidate_tax}))
			{
			print "\twarning, taxon($candidate_tax) was come across on previous node. ignoring.\n";
			}else{
			#print "child $i of node $next has tax label:$candidate_tax\n";
			#print "\t\tthis taxon has new IDs to be added\n";

			# these are new IDs, comma seperated. will have comma terminating string.
			my $new_members = $assign_these_new_terminals_to_taxon{$candidate_tax};# . ",";
			$new_IDs_going_to_this_node .= $new_members;
			$new_members =~ s/\,$//;
			my @new_members_to_node = split /\,/ , $new_members;
			print "\tcount new members to this node:" , scalar @new_members_to_node , "\n";
			};

		$assigned_already{$candidate_tax} = 1;
		}else{#if(exists($assign_these_new_terminals_to_taxon{$candidate_tax}))
		print "\tNO new members to assign to this name\n";
		}

	};#foreach my $candidate_tax(@split_node_tax_labels)






# constrinat is made of the new taxa, and the list of terminasl in backbone tree, so get the terminals:
$all_terminals_from_this_node = "";
# does this:$all_terminals_from_this_node .= "$test\t";

#######################################
get_terminals_from_this_node2($next);#
#######################################





unless($next eq $root_node)
	{



# 3 components to each constraint: 1)the taxa monophyletic 2) the taxa excluded 3) the taxa that float.
# err on the side of floating!

my %floaters = ();
my %constrained = ();

# all backbone IDs:		%terminals
# backbone IDs to constrain:	$all_terminals_from_this_node


my @all_backbone_IDs = keys %terminals;@all_backbone_IDs = sort @all_backbone_IDs;
foreach my $backbone_ID(@all_backbone_IDs)
	{
	#print "backbone_ID:$backbone_ID\n";	
	if($all_terminals_from_this_node =~ /$backbone_ID/)
		{
		#print "\tconstrain MP\n";
		$constrained{$backbone_ID}=1;# taxa in the backbone tree decended from this node will be constrained MP
		}else{
		$constrained{$backbone_ID}=2;# taxa in backbone tree not decended from this node will be excluded from MP
		#print "\texclude from MP\n";#die;
		};
	};

#$all_terminals_from_this_node

my @all_new_IDs = keys %query_IDs;
foreach my $new_ID(@all_new_IDs)
	{
	my $check_overlapping = $constrained{$new_ID};
	#print "new_ID:$new_ID check_overlapping:$check_overlapping\n";
	unless($check_overlapping =~ /\d/)
	{
	if($new_IDs_going_to_this_node =~ /$new_ID/)
		{
		#print "\tconstrain MP\n";
		$constrained{$new_ID}=1;# new taxa assigned to this node will be constrained MP
		}else{
		#print "\tfloat\n";
		$constrained{$new_ID}=3;# new taxa not assigned to this node will float
		};
	}};



my @everyone = keys %constrained;@everyone = sort @everyone;

my $constrain_string="";
my $float_string="";
my %check_assigned = ();
$constraint_number++;


foreach my $member(@everyone)
	{
	#print "member:$member\n";
	my $contr = $constrained{$member};
	my $fasta_index;
	if($index_of_fastaID{$member} =~ /\d/)
		{
		$fasta_index = $index_of_fastaID{$member};
		if(exists($check_assigned{$fasta_index}))
			{die "\nerror 2218. member:$member fasta_index:$fasta_index\n"};
		$check_assigned{$fasta_index}=1;
		}else{die "\nerror 2190. constr:$contr\tmem:$member\t\n"}
	#print "$contr\t$member\tfile index:$fasta_index\n";

	if($contr == 1){$phylip_constraints{$member} .= "1";$constrain_string .= "$fasta_index "}
	if($contr == 2){$phylip_constraints{$member} .= "0"};
	if($contr == 3){$phylip_constraints{$member} .= "-";$float_string .= "$fasta_index "}
	}

$constrain_string	=~ s/\s$//;
$float_string		=~ s/\s$//;

my $tnt_command;
$tnt_command = "[ $constrain_string ($float_string) ]";

$tnt_command_string .= "force + $tnt_command, ";
#print TNT_FILE_HANDLE "$tnt_command_string";


#		$index_of_fastaID{$id} = $count_diet_IDs;
# tnt format, 
#	# constraint where 2 taxa 'float'
#	force + [ a b c  (d e) ]
#tnt64 mxram 2000, rseed 1,p default-Second.nex timeout 0:20:00, echo=, taxname=, force + [17 18 19 20], constrain =,
#
#$tnt_command_string
	# crashes unless you have single quotes around force command, due to parentheses confusing bash.
#	tnt64 mxram 2000, rseed 1,p tnt_input.nex timeout 0:00:30, echo=, taxname=, 'force + [ 14 46 67 (1 2) ],' constrain =, log tnt_logfile,  hold 2000, xmult=replication 100 drift 10 ratchet 10 hold 10, best, majority *, tchoose / , export *tnt_output, log / ,quit;
#print "backbone constraint (list of all_terminals_from_this_node):$all_terminals_from_this_node\n";
#print "new_IDs_going_to_this_node:$new_IDs_going_to_this_node\n";


# fasttree constraint file.
# 5 2
#A        00
#B        0-
#C        10
#D        11
#E        -1




	};#unless($next eq $root_node)




for $i(0 .. $child_nodes)# get child nodes as normal
	{
	push @next1 , $nodes{$next}{$i};
	#print "child node:$i node ID:$nodes{$next}{$i}\n";
	};

for my $index(0 .. $#next1)
	{
	my $test = $next1[$index];#print "test:$test\n";

	if($test =~ /^INTERNAL_NODE_/)
		{
		#####################################
		get_terminals_from_this_node31($test  );#	# recurse
		#####################################
		}
	}

return();

}


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





sub traverse_backbone_tree_and_collapse
{

my $next = shift;

# this has species names, species groups, genera, families etc. and has the outgroup on first call.
my $all_taxnames_for_terminals = $all_taxnames_for_terminals_decended_from_this_node{$next};
#print "all_taxnames_for_terminals:($all_taxnames_for_terminals)\n";


# taxonomic name assigned to node of backbone tree:
my $label = $node_label{$next};# just one name in this
# or this one:
#my $label = $nodes{$next}{node_taxonomic_label};
#print "\nNEW NODE ($next), label:$label\n";


my @non_terminal_child_nodes = ();my @child_nodes_test = ();
my $all_terminals_from_both = "";my @taxon_for_child_nodes = ();
my $lowest_taxname_for_node = "";my $lineage_for_node = "";
$all_terminals_from_this_node ="";


#my @tax_labels_current_node_of_backbone_tree = split /\./ , $label;



#foreach my $taxlab(@tax_labels_current_node_of_backbone_tree)
#	{
my $taxlab = $label;
my $tax_number = $ncbi_tax_number_for_this_species{$taxlab};

# all taxa at rank below the name:
my $child_nodes = $ncbi_nodes{$tax_number}{child_nodes};$child_nodes =~ s/^\t//;

# $ncbi_nodes{$tax_number}{child_nodes};
# $ncbi_nodes{$child}{name};

#print "
#label assigned to current node of backbone:$taxlab, NCBI tax_number:$tax_number
#tax_number:$tax_number has child_nodes:$child_nodes\n";

my @child_nodes_array = split(/\t/, $child_nodes);
my $current_constraint = "";my $current_constraintB = "";

my $mpA =0;my $mpB =0;my $data_missing_reached=0;
my $count_floaters = 0;my $count_constrained = 0;

# each child node (taxonomic name) derived from name that was given to current node of tree
foreach my $child(@child_nodes_array)
	{
	my $name_string = $ncbi_nodes{$child}{name};$name_string =~ s/\s+/_/g;	
	# print "\tchild taxon ID:$child name:$name_string\n";

	# for each child taxon (1 rank below) derived from the 
	# taxonomic label assigned to the current node of the backbone tree.

	my $found_in_backbone_decendent=0;
	if($all_taxnames_for_terminals =~ /\t$name_string\t/)
		{
		# child taxon is found in the lineage of terminals decended from this backbone node			
		# so a monophyletic constraint may be made			
		$found_in_backbone_decendent =1;#print "\tIS in list of taxnames of decendents\n";

		#print "\t\t\ttax found in a terminal of this node of backbone tree\n"
		}else{
		#print "\tNOT in list of taxnames of decendents\n";
		#print "\t\t\ttax NOT found in a terminal of this node of backbone tree\n";
		};

	# how many times has this taxon been assigned, ANYWHERE on the backbone tree
	my $number_anywhere_in_backbone = $number_of_nodes_to_which_this_taxon_has_been_assigned{$name_string};

	# is this taxon found in the fasta file, if not, ignore.
	my $how_many_in_fasta = $count_tax_in_new_data{$name_string};

	# if child taxon if found in backbone tree and new fasta members, but is not in decendents of current node, 
	# will constrained outside of monophyly
	if($number_anywhere_in_backbone >= 1 && $how_many_in_fasta >= 1 && $found_in_backbone_decendent == 0)
		{$current_constraint .= "$name_string:0\n";$mpA=1;
		$count_constrained += $how_many_in_fasta;
	#	print "\t\tmpA:$mpA, constrain $child:$name_string ($how_many_in_fasta) outside MP\n"
		}

	# constrained MP
	if($number_anywhere_in_backbone >= 1 && $how_many_in_fasta >= 1 && $found_in_backbone_decendent == 1)
		{$current_constraint .= "$name_string:1\n";$mpB=1;$current_constraintB .= "$name_string-";
		$count_constrained += $how_many_in_fasta;
	#	print "\t\tmpB:$mpB, constrain $child:$name_string ($how_many_in_fasta) inside MP\n"
		}

	# no information on where these should be placed (not on backbone). so they float where they want.
	if($number_anywhere_in_backbone <1)
		{
		if( $how_many_in_fasta >= 1)
			{$current_constraint .= "$name_string:-\n";
			$data_missing_reached =1;
			$count_floaters += $how_many_in_fasta;
		#	print "\t\tnot in backbone, BUT IS in fasta ($how_many_in_fasta), unknown for constraining\n";
			}else{
		#	print "\t\tnot in backbone, AND NOT in fasta, so, inconsequential.\n";
			}
		};
#	print "\t$child:$name_string how_many_in_fasta:$how_many_in_fasta\n";
#	print "\tnbr_anywhere_in_backbone:$number_anywhere_in_backbone found_in_backbone_decendent:$found_in_backbone_decendent\n";


	};#foreach my $child(@child_nodes_array)


my $constraint_length = length($current_constraint);


# mpA and mpB check that there are taxa defined within and external , i.e. a monophyly to be made.

#####################################################################

my $test1 = 0; if($mpA ==1 && $mpB == 1 && $count_floaters == 0){$test1 = 1};

my $current_constraint_string = "";
if($node_monophylys{$next} =~ /\w/ || $test1 == 1)
	{
	print "\ndefine cosntraint for BACKBONE NODE NUMBER:$next\n";	
	};

my $monophyletic_taxon_defined_from_this_node = "NA";
if($node_monophylys{$next} =~ /\w/)
	{
	$monophyletic_taxon_defined_from_this_node = $node_monophylys{$next};
	print "\tmonophyletic_taxon:$monophyletic_taxon_defined_from_this_node\n";
	$current_constraint_string = $monophyletic_taxon_defined_from_this_node;
	};

 if($mpA ==1 && $mpB == 1 && $count_floaters == 0)
	{
	$current_constraintB =~ s/-$//;
	print "\trelational constraint: $current_constraintB\n";
	$current_constraint_string = $current_constraintB;
	$ft_constraints{$current_constraint}=1;	# seems not used
	};

# some duplicates:
if($constraints_written{$current_constraint_string} == 1) # %constraints_written does seem to be used.
	{$current_constraint_string = ""};
$constraints_written{$current_constraint_string} = 1;

#####################################################################



my $child_nodes = $child_counts{$next};
my @next1 = ();

for $i(0 .. $child_nodes)
	{
	# specify new array of the child nodes as normal for recursing
	my $push_node = $nodes{$next}{$i};
	push @next1 , $push_node;	
	}


my $swap_string 	= "";
$join_the_child_nodes 	= join ',', @next1;

if($current_constraint_string =~ /\w/)
	{	
	$swap_string 		= "($join_the_child_nodes)$current_constraint_string";
	
	}else{
	$swap_string 		= "$join_the_child_nodes";
	
	};

$yet_another_newick_string 	=~ s/$next(\W)/$swap_string$1/;


# alternative newick string that is readable at this stage:

if(exists($child_counts{$next}))
	{
	#$new_newick_string2 	=~ s/$next(\W)/$swap_string$1/;
	}else{
	# instead of replacing node (which is now a species name) with child nodes (which there arnt any)
	# replace node with lineage
	#$new_newick_string2 	=~ s/$next(\W)/$label_string_to_assign$1/;
	}


if(exists($child_counts{$next}))
	{
	for my $index(0 .. $#next1)
		{
		my $test = $next1[$index];

		#####################################
		traverse_backbone_tree_and_collapse($test );#	# recurse
		#####################################
		}
	};

return();

}






#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################





sub get_all_taxa
{


my $tax_list = shift;
my @tax_array = split /\t/ , $tax_list;
my $all_tax_names_for_these_terminals = "";

unless($tax_list =~ /./){die "\nerror 2830, why has sub get_all_taxa been called with nothing\n"};

foreach my $tax(@tax_array)
	{
	my $test_taxonomy = $complete_lineage_for_this_species{$tax};
	my $genus_name = $tax;$genus_name =~ s/^([A-Z][a-z]+)_.+/$1/;
	unless($test_taxonomy =~ /\w/){$test_taxonomy = $complete_lineage_for_this_species{$genus_name}};

	#print "\tspecies:$tax complete taxonomy:$test_taxonomy2\n";

	while($test_taxonomy =~ s/^([^:]+):(\w+)//)
		{
		my $current_rank = $1;my $current_taxname = $2;
		unless($all_tax_names_for_these_terminals =~ /\t$current_taxname\t/)
			{$all_tax_names_for_these_terminals .= "	$current_taxname	";
			}
		}
	}
#print "\ncount tax_array:" , scalar @tax_array , ", all_members_have_this_name:$all_members_have_this_name\n";


return($all_tax_names_for_these_terminals);


}






#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


sub test_for_monophylies
{
my $tax_list = shift; # input is list of all terminals from current node
my @tax_array = split /\t/ , $tax_list;
my %all_tax_names_for_these_terminals = ();
my $number_input = scalar @tax_array;
print "\nsub test_for_monophylies, number terminals input:$number_input\n";


#count_taxa_represented_in_terminals


foreach my $tax(@tax_array)
	{
	my $test_taxonomy = $complete_lineage_for_this_species{$tax};
	unless($test_taxonomy =~ /\w/)
		{
		my $genusname = $tax;$genusname =~ s/^([A-Z][a-z]+)_.+/$1/;
		$test_taxonomy = $complete_lineage_for_this_species{$genusname};
		
		if($test_taxonomy =~ /\w+/)
			{
			#print "\tphew, found using genus name ($genusname), continuing ...\n";
			}else{
			die "\nerror 2889, no taxonomic informartion found for ($tax) nor ($genusname)\n";		
			}
		};
	#print "\tspecies:$tax\n";# complete taxonomy:$test_taxonomy\n";
	while($test_taxonomy =~ s/^([^:]+):(\w+)//)
		{
		my $current_rank = $1;my $current_taxname = $2;
		$all_tax_names_for_these_terminals {$current_taxname}++;
	#	print "\t\tcurrent_rank:$current_taxname:" , $all_tax_names_for_these_terminals{$current_taxname}, "\n";
		}
	}



#count_taxa_represented_in_terminals

# for a given node of the backbone tree, there may be several tax names for which the terminlas are monophyletic,
# only the lowest (tip-ward) rank is valid,
# for example just 2 terminals are input into the monophyly test sub,
# these 2 belong to Chrysidini tribe, and the only representatives of the superfamily Chrysidoidea,
# thus each of Chrysidini, Chrysidinae, Chrysididae, Chrysidoidea
# would 'monophyletic' in loosest test, but for Chrysidinae and above this would be too permissive a test,
# as for the most part there would be only one (Chrysidini in this case) of 
# presumably several fundamental representatives (other tribes),
# so, here need to get rank of each, and select lowest where several are present.

my @all_taxnames_from_terminals_of_current_node = keys %all_tax_names_for_these_terminals;
@all_taxnames_from_terminals_of_current_node = sort @all_taxnames_from_terminals_of_current_node;
my $current_no_taxnames = scalar @all_taxnames_from_terminals_of_current_node;
my $monophyly_found=0;
my $putative_monophhyly_rank = 0;my $putative_monophhyly_name = "NA";

foreach my $tax(@all_taxnames_from_terminals_of_current_node)
	{
	my $count_total_in_backbone_terminals = $count_taxa_represented_in_terminals{$tax};
	my $count_decended_from_current_node = $all_tax_names_for_these_terminals{$tax};
	my $how_many_terminals_from_current_node = scalar @tax_array;
	my $rank_number = $rank_codes{$tax};
	unless($rank_number =~ /\d/)
		{$rank_number = 1;print "wanriing, couldtn find rank number for tax:$tax\n"};

#	print "tax:$tax rank_number:$rank_number total_in_backbone_terminals:$count_total_in_backbone_terminals\n";

	if($count_total_in_backbone_terminals == $count_decended_from_current_node && 
		$count_total_in_backbone_terminals == $how_many_terminals_from_current_node)
		{
		# this was too liberal a definition:
	#	$monophylys{$tax}=1; 

		$monophyly_found = 1; print "putative monophyly:$tax\n";

		if($rank_number >= $putative_monophhyly_rank)
			{
			$putative_monophhyly_rank = $rank_number;
			$putative_monophhyly_name = $tax;
			};
		
		}
	}

 if($monophyly_found == 1 )
	{
	$monophylys{$putative_monophhyly_name} = 1;
	print "DEFINED monophyly:$putative_monophhyly_name\n";
	};


# preiovusly was this, though seems the returned varable not used
# return($all_tax_names_for_these_terminals);



# this is now required retunred:
return($putative_monophhyly_name);







};





#####################################################################################################
#
#
#
#####################################################################################################


sub place_specieslevel_data_into_taxon_constraints
{

my $newick_constraint_string 	= shift;
my $newick_length 		= length($newick_constraint_string);

print "\n
sub place_specieslevel_data_into_taxon_constraints
newick_length:$newick_length
\n";

open(VERBOSE_LOG, ">backbone_constraints_newick_verbose_LOG");


my $interal_node	= 0;

while ($newick_constraint_string =~ s/\(([^\(\)]+)\)([a-zA-Z\-]*)/INTERNAL_NODE_$interal_node/) # processed sumtrees output from a load of raxml boots
	{

	my $node = $1;my $internal_label = $2; my $nodeID = "INTERNAL_NODE_$interal_node"; #print "nodeID:$nodeID node:$node\n";
	my $node_length = length($node);
	# contains multi-furcations.
	my @child_nodes = split /\,/ , $node;my $number_child_nodes = scalar @child_nodes;

	

	push @child_nodes , $internal_label;
	$all_internal_node_strings_assigned_to_input_constraint{$internal_label} = 1;


#	print "\nfirst adding new members of child nodes\n";
	for $i ( 0 .. $#child_nodes )
		{
		print VERBOSE_LOG "nodeID:$nodeID node_length:$node_length count_child_nodes:$#child_nodes current_child_no:$i, child:$current_child, count added $count_added_to\n";

		#print "$child_nodes[$i]\n";
		my $current_child = "-$child_nodes[$i]-";
		if($current_child =~ /[\,]/){die "\nerror 3084:$current_child\n"};
		# $new_member_assigned_to_constraint_taxon { $new } = $current_member_to_be_assigned;
		# my @testkeys = keys %new_member_assigned_to_constraint_taxon;
		my @add = ();

		foreach my $new ( @testkeys )
			{
			my $is_assigned_to = $new_member_assigned_to_constraint_taxon { $new };
			print VERBOSE_LOG "\tnew:$new is_assigned_to:$is_assigned_to\n";
			if($current_child =~ /-$is_assigned_to-/)
				{
				unless($this_added_already{$new} == 1) 
					{
					print VERBOSE_LOG "\t\tADDING\n";
					 push @add , $new 
					};
				$this_added_already{$new} = 1;
				};
			};
		my $count_added_to = scalar @add;
		$current_child =~ s/^-(.+)-$/$1/;
		if($count_added_to >= 1){$add_these_to_this{$current_child} = join ',' , @add};

		};

#	print "node:$interal_node, number_child_nodes:$number_child_nodes int_label:$internal_label\n\tchild nodes:@child_nodes\n";

	# node:4, int_label:Thysanoptera
	# 	child nodes:Thripinae
	# node:5, int_label:Psocoptera-Thysanoptera-Phthiraptera
	# 	child nodes:Psocoptera-Phthiraptera INTERNAL_NODE_4





	$interal_node++;
	}


print "
remaining \$newick_constraint_string (should be nothing):$newick_constraint_string
";




# finally, string replace newick contraints with new data

my @keys = keys %add_these_to_this;

# open(LOG2 , ">$outfile_prefix.bcn.LOG") || die "\nerror 3441\n";

foreach my $key(@keys)
	{
	my $replace_with = $add_these_to_this{$key};
	my $count_members = scalar split /\,/ , $replace_with;

	my $screen_print = $replace_with;$screen_print =~ s/^(.{120}).+/$1........../;
	#if(length($replace_with)<= 1000){
		print "replcace assign tax $key with $count_members members:$screen_print\n";
	#	print LOG2 "replcace assign tax $key with $count_members members:$screen_print\n";
print VERBOSE_LOG "replcace assign tax $key with $count_members members:$screen_print\n";

	#};

		# put comma after bracket otherwise will end up as node label:
	if($yet_another_newick_string =~ s/([\)])$key([\(\)\,])/$1,$replace_with$2/) # bugfix 20160928
		{}elsif($yet_another_newick_string =~ s/([\(\)\,])$key([\(\)\,])/$1$replace_with$2/)
		{
		my $test1 = $1;my $test2 = $2;#print "\tpreceednig char:$test1 following char:$test2\n";
		}else{
		print "\nerror 3441???\n";
		};

	};

# close LOG2;
close VERBOSE_LOG;

# correct clade with onyl one tip. a few node lables will be removed by this, but most remain.

# (Challia)Orthopteroidea-Neoptera_incertae_sedis,
$yet_another_newick_string =~ s/\(([A-Z][a-z]+)\)[a-zA-Z\-\_]+/$1/g;
# ERROR: Expecting ',' in tree; found: character ')'
# chnura)),((((Challia),(Mesocapnia,Aptero
$yet_another_newick_string =~ s/\(([a-z]+)\)\,/$1,/ig;






} # sub

#####################################################################################################
#
#
#
#####################################################################################################



sub backbone_constraints
{


	# THIS SUB NOT RUN



# this didnt work, if backbone trees were perfect then maybe it would

print "
sub backbone_constraints
";

my @all_backbone_terminals = keys %terminals;@all_backbone_terminals = sort @all_backbone_terminals;


	# THIS SUB NOT RUN


foreach my $termnial(@all_backbone_terminals)
	{
	my $parent = $nodes{$termnial}{parent};
#	print "termnial:$termnial parent:$parent\n";
	my $sister_node;
	for my $i(0,1)
		{
		my $child_node_of_parent = $nodes{$parent}{$i};
		if($child_node_of_parent eq $termnial)
			{}else{
			$sister_node = $child_node_of_parent;
			};
	#	print "\ti:$i child_node_of_parent:$child_node_of_parent\n";
		};
	unless($termnial =~ /\w/ && $parent =~ /\w/ && $sister_node=~ /\w/){ die "\nerror 3552\n"};
	print "\tsister_node:$sister_node\n";

	my $most_inclusive_taxa_of_sister_node = "";
	if($sister_node =~ /INTERNAL_NODE_/)
		{
		$all_terminals_from_this_node = "";
		get_terminals_from_this_node2($sister_node, 0);#
		$most_inclusive_taxa_of_sister_node = get_shared_taxa($all_terminals_from_this_node);#
		unless($most_inclusive_taxa_of_sister_node =~ s/.+\t(\w+)\t$/$1/)
			{
			print "\nwarning, expecting lineage here, algorthim won work for this terminal\n";
			$most_inclusive_taxa_of_sister_node = "NA";
			};
		}else{
		$most_inclusive_taxa_of_sister_node = $sister_node;
		};
	print "\tmost_inclusive_taxa_of_sister_node:$most_inclusive_taxa_of_sister_node\n";

	# THIS SUB NOT RUN


	################################################################################

	unless($most_inclusive_taxa_of_sister_node eq "NA")
		{
		# for genus of current terminal, walk up the taxonomic ranks, 
		# until taxa is shared by sister of backbone tree.
		
		my $test_species = $termnial;
		my $test_taxonomy = $complete_lineage_for_this_species{$test_species};
		unless($test_taxonomy =~ /\w+/)
			{
			my $genusname = $test_species;$genusname =~ s/^([A-Z][a-z]+)_.+/$1/;
			$test_taxonomy = $complete_lineage_for_this_species{$genusname};
			if($test_taxonomy =~ /\w+/)
				{}else{
				if($test_species =~ /optera_[A-Z]/)
					{print "\nlooks like user format was not processed to an established tax name\n"};
				die "\nerror 3781, no taxonomy found for species ($test_species) nor genus ($genusname). test_species:($test_species) tax_list:($tax_list) \nquitting ...\n";
				};
			};


	# THIS SUB NOT RUN


		print "test_taxonomy:($test_taxonomy)\n";
		my $upper_limit=0;my $upper_taxon = "NA";
		while($test_taxonomy =~ s/\s+(\w+):(\w+)\s+$/ /)
			{
			my $current_rank = $1;my $current_taxname = $2;

			unless($upper_limit == 1)
			{
			print "\t\tNEW rank:$current_rank taxname:$current_taxname\n";
			my $ncbi_taxnumber = $ncbi_taxnumber_for_taxname{$current_taxname};
			unless($ncbi_taxnumber =~ /[\d\w]/){die "\ndidnt find tax number for name\n"};

			# traversing heirachy from this node, 
			# checking no offspring are sister taxa
			$all_taxa_from_this_node="";
			traverse_heirachy($ncbi_taxnumber);
		#	print "all_taxa_from_this_node:$all_taxa_from_this_node\n";
			if($all_taxa_from_this_node =~ /\t$most_inclusive_taxa_of_sister_node\t/)
				{
				$upper_limit = 1;
				};
			if($upper_limit == 0){$upper_taxon = $current_taxname};

			};
			};
			print "highest taxon not shared by backbone sister:$upper_taxon\n";
		$backbone_terminal_rename{ $termnial } = $upper_taxon;
		};
};



	# THIS SUB NOT RUN


foreach my $termnial(@all_backbone_terminals)
	{
	my $replace_backbone_genus_with_string = "";

	if($backbone_terminal_rename{ $termnial } =~ /\w/)
		{
		$replace_backbone_genus_with_string = $backbone_terminal_rename{ $termnial };

		if($process_backbone_tree_terminal_IDs == 0)
			{
			$second_backbone_copy_newick =~ s/$termnial([^a-zA-Z])/$replace_backbone_genus_with_string-$termnial$1/;
			}else{
			$second_backbone_copy_newick =~ s/$termnial([\_])/$replace_backbone_genus_with_string-$termnial$1/;
			};		

		};


	};

	# THIS SUB NOT RUN

close TEST_FILE;


open (OUT6999 , ">$outfile_prefix.less_basic_constraint_tree") || die "\n\n";
print OUT6999 "$second_backbone_copy_newick\n";
close OUT6999;



	# THIS SUB NOT RUN




die "";




};



#####################################################################################################
#
#
#
#####################################################################################################






sub traverse_heirachy
{
my $current_node = $_[0]; #my $current_node_taxstring = $_[1];# $current_node_taxstring to be deprecated
my $child_nodes = $ncbi_nodes{$current_node}{child_nodes};$child_nodes =~ s/^\t//;
my @child_nodes_array = split(/\t/, $child_nodes);

if($current_node == $starting_node)
	{
	my $rank = $ncbi_nodes{$starting_node}{rank};$rank =~ s/\s/_/g;
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
	my $rankcode = $ncbi_nodes{$child}{rank_code};$rank_codes{$name_string} = $rankcode;

	my $child_complete_lineage = $ncbi_nodes{$current_node}{complete_lineage} . "$rank:$name_string ";#prev assigned $name
	my $originalname = $ncbi_nodes{$child}{name};$originalname =~ s/[\s\t]/_/g;
	my $name_assignment_to_taxnumber = "";



	$all_taxa_from_this_node .= "\t$originalname\t"; ######################    <-   this is the only thing happening here



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

		###########################################
		traverse_heirachy($child );#
		###########################################
	}}

	
}#sub traverse_heirachy





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



sub backbone_constraints2
{


# attempt 2 .. 'BOTTOM UP' constraining. easily maintains integrity of the backbone tree but may limit number of barcodes which can be assigned



#
#
#
#
# taxonomic assignment in backbone tree, pseudocode:
# for each backbone terminal
#	retreive its taxonomic lineage
#	for each parent taxon in lineage
#		if taxon observed in a different terminal of the phylogeny, store lower taxon name with terminal name, break. 
#
# store names of backbone ids with no barcodes for pruning
#
#
#

print "
sub backbone_constraints2\n
";

my @all_backbone_terminals = keys %terminals;@all_backbone_terminals = sort @all_backbone_terminals;

print BCN_LOG "\nrunning algorithm to assign higher taxa to backbone terminals ...\n";

foreach my $termnial(@all_backbone_terminals)
	{
	my $new_lineage_string = "";
	$terminal_count++;
	$terminal_print2++;if($terminal_print2 <= 6){print "$terminal_count of $#all_backbone_terminals\n"}
	if($terminal_count =~ /00$/)
		{$terminal_print++;if($terminal_print >= 7 && $terminal_print <= 12){print "$terminal_count of $#all_backbone_terminals\n"}
		};

		my $test_species = $termnial;my $test_taxonomy = $complete_lineage_for_this_species{$test_species};
		unless($test_taxonomy =~ /\w+/)
			{
			my $genusname = $test_species;$genusname =~ s/^([A-Z][a-z]+)_.+/$1/;$test_taxonomy = $complete_lineage_for_this_species{$genusname};
			if($test_taxonomy =~ /\w+/)
				{}else{
				if($test_species =~ /optera_[A-Z]/){print "\nlooks like user format was not processed to an established tax name\n"};

				print "\nerror 4046, no taxonomy found for species ($test_species) nor genus ($genusname). test_species:($test_species) tax_list:($tax_list) \n";
				print BCN_LOG "warning, no taxonomy found for terminal $test_species, nor genus ($genusname).\n";
				if($genusname =~ /^Beetles$|^Flies$|^Moths$|^Wasps$/){}else{
					#die ""
					};

				};
			}; 
	#	 print "TERMINAL:$termnial test_taxonomy:($test_taxonomy)\n";

		my $upper_limit=0;my $upper_taxon = "NA";my $break_loop = 0;


# HACK, dont replace outgroup with a load of its relatives
#	if($termnial =~ /Boreus_hyemalis/){$upper_taxon = "Boreus_hyemalis";$new_lineage_string .= "$current_taxname;";$break_loop = 1};
#	if($termnial =~ /Stylops_melittae/){$upper_taxon = "Stylops_melittae";$new_lineage_string .= "$current_taxname;";$break_loop = 1};
#	if($termnial =~ /Corydalus_cornutus/){$upper_taxon = "Corydalus_cornutus";$new_lineage_string .= "$current_taxname;";$break_loop = 1};

#	if($termnial =~ /Macropis_europaea/){$upper_taxon = "Macropis_europaea";$new_lineage_string .= "$current_taxname;";$break_loop = 1};

	if($retain_terminal =~ /\w/)
		{
		if($termnial =~ /$retain_terminal/){$upper_taxon = $retain_terminal;$new_lineage_string .= "$current_taxname;";$break_loop = 1};
		};




		while($test_taxonomy =~ s/\s+(\w+):(\w+)\s+$/ /)
			{
			my $current_rank = $1;my $current_taxname = $2;
		#	print "\trank:$current_rank tax:$current_taxname\n";

		#	if($termnial =~ /Micropterix/){print "\trank:$current_rank tax:$current_taxname\n"};

			unless($break_loop == 1)
			{
			foreach my $termnial2(@all_backbone_terminals)
				{
				unless($termnial2 eq $termnial)
					{
					my $test_species2 = $termnial2;
					my $test_taxonomy2 = $complete_lineage_for_this_species{$test_species2};
					unless($test_taxonomy2 =~ /\w+/)
						{
						my $genusname2 = $test_species2;$genusname2 =~ s/^([A-Z][a-z]+)_.+/$1/;
						$test_taxonomy2 = $complete_lineage_for_this_species{$genusname2};
						if($test_taxonomy2 =~ /\w+/)
							{}else{
							if($test_species2 =~ /optera_[A-Z]/){print "\nlooks like user format was not processed to an established tax name\n"};
							$backbone_terminals_for_which_taxonomy_not_found{$test_species2} =1;

							$error_printing++;
							if($error_printing < 5){
							print "\nerror 4078, terminal $terminal_count. no taxonomy found for $test_species2 nor genus ($genusname2). test_species:($test_species) tax_list:($tax_list) ...\n";
									}elsif($error_printing == 5){print "\nnot printing more errors\n"};

							if($genusname2 =~ /^Beetles$|^Flies$|^Moths$|^Wasps$/)
								{
								}else{
								# die ""
								};
							};
						};#print "\ttermnial2:$termnial2 test_taxonomy2:($test_taxonomy2)\n";

					#######################
					while($test_taxonomy2 =~ s/\s+(\w+):(\w+)\s+$/ /)
						{
						my $current_rank2 = $1;my $current_taxname2 = $2;
						if($termnial =~ /Micropterix/ && $test_taxonomy2 =~ /Lepidopters/)
							{
						#	print "\tcurrent_rank2:$current_rank2 current_taxname2:$current_taxname2\n"
							};
	
						if($current_taxname eq $current_taxname2)
							{
							$break_loop = 1;
							};

						};
					#####################

					};
				};# foreach my $termnial2(@all_backbone_terminals)

			unless($break_loop == 1)
				{
				$upper_taxon = $current_taxname;
				$new_lineage_string .= "$current_taxname;";
			#	if($termnial =~ /Micropterix/){print "\tupper_taxon:$upper_taxon\n"};

				};

			};# unless($break_loop == 1) # this stops at upper taxon
			}; # while($test_taxonomy =~


	

		#	print "terminal:$termnial new_taxon:$upper_taxon\n";
		$backbone_terminal_rename{ $termnial } = $upper_taxon;
		$backbone_terminal_rename2{ $termnial } = $new_lineage_string;

			
	}; # foreach my $termnial(@all_backbone_terminals)


print "done ... looking through terminals again, replacing strings ... \n";

print BCN_LOG "\n";


##############################################


# this bit just stores terminal name with replcing taxon, and couple other things. removes terminals from list in which higher taxon not inferred.
# ignore for pseudocode


open(LOG_OUT , ">Backbone_Constraints_Newick.list_of_constraints") || die "";

$terminal_count = 0;$terminal_print2 = 0;$terminal_print = 0;
foreach my $termnial(@all_backbone_terminals)
	{
	$terminal_count++;
	$terminal_print2++;
	if($terminal_count =~ /00$/){$terminal_print++;if($terminal_print >= 7 && $terminal_print <= 12){print "$terminal_count of $#all_backbone_terminals\n"}};

	my $replace_backbone_genus_with_string = "";

	if($backbone_terminal_rename{ $termnial } =~ /\w/)
		{

		$replace_backbone_genus_with_string = $backbone_terminal_rename{ $termnial };
		$backbone_substrings_original_length += length($termnial);
		$backbone_substrings_replace_length += length($replace_backbone_genus_with_string);

		if( $process_backbone_tree_terminal_IDs == 0 )
			{
			unless($second_backbone_copy_newick =~ s/([\,\(])$termnial([^a-zA-Z])/$1$replace_backbone_genus_with_string-$termnial$2/)
				{
				die "\nError 4198, could not replace string in newick:$termnial .. quitting.\n"
				};
			
			}else{
			unless($second_backbone_copy_newick =~ s/([\,\(])$termnial([\_])/$1$replace_backbone_genus_with_string-$termnial$2/)
				{
				die "\nError 4204, could not replace string in newick:$termnial .. quitting.\n"
				};
			};


		unless($backbone_terminal_rename{ $termnial } eq "NA")
			{
		#	$store_all_constraints{$replace_backbone_genus_with_string} = "$replace_backbone_genus_with_string-$termnial"
		# try and make simple as possible , replace with new tax lists, on original newick string, not taxon assigned string:
			$store_all_constraints{$termnial} = $replace_backbone_genus_with_string;
			};

		if($terminal_print2 <= 6)
			{
			print "$terminal_count of $#all_backbone_terminals:$termnial\n";
			print "  replace terminal ($termnial) with higher taxon ($replace_backbone_genus_with_string)\n";
			};



		}; # 	if($backbone_terminal_rename{ $termnial } =~ /\w/)

	 print LOG_OUT "$termnial\t$replace_backbone_genus_with_string\n";

	};

close TEST_FILE;
close LOG_OUT;

print "
backbone_substrings_original_length:$backbone_substrings_original_length
backbone_substrings_replace_length:$backbone_substrings_replace_length
";

#############################################################################




#
# pseudocode
# for each backbone terminal in which higher taxon inferred,
#	retreive list of barcode IDs belonging to taxon,
#	replace posiiton in newick string with comma-seperated list
#
#


print "
replacing higher taxon names with barcode IDs 
";


open (OUT69999 , ">$outfile_prefix.less_basic_constraint_tree.backbone_taxa_assigned") || die "\n\n";
print OUT69999 "$second_backbone_copy_newick\n";
close OUT69999;

# open(BCN_LOG, ">>backbone_constraints_newick_LOG") || die "";


if($graft_source_trees_list =~ /\w/)
	{
	open(GRAFT_LOG, ">>grafting_LOG.txt") || die "";
	print GRAFT_LOG "\ninferring higher taxa for backbone terminals\n";
	};

my @list_constraints = keys %store_all_constraints; @list_constraints = sort @list_constraints;

my $cosntr_number = 0;
$grafts_made =0; $grafted_count=0;

foreach my $termnial(@list_constraints) # this list of terminals is not neccessary complete,
	{					# will omit those in which upper_taxon inference did not work

	my $constreaint = $store_all_constraints{$termnial}; # hash key is backbone terminal, entry is upper taxon assigned to the terminal


	$cosntr_number++; 	# print "\nNEW terminal:$termnial constraint:$constreaint new members($fasta_ids_for_each_taxon{$constreaint})\n";
	if( $fasta_ids_for_each_taxon{$constreaint} =~ /\S/ )	# one or more
		{
		my $new_members_included_in_constraint = $fasta_ids_for_each_taxon{$constreaint};
		$new_members_included_in_constraint =~ s/^\t+//;$new_members_included_in_constraint =~ s/\t+$//;
		$new_members_included_in_constraint =~ s/\t\t+/	/g;
		my @split_memers = split /\t+/ , $new_members_included_in_constraint;
		my $count_new_members = scalar @split_memers;
		print BCN_LOG "terminal count $cosntr_number, termnial:$termnial, higher taxon_assigned:$constreaint, COI matching:$count_new_members members:$new_members_included_in_constraint\n";
		print GRAFT_LOG "\nCONSTR number:$cosntr_number, termnial:$termnial, taxon_assigned:$constreaint, count_new:$count_new_members new_members_incl:$new_members_included_in_constraint\n";

#		my $find_string = $store_all_constraints{$constreaint};
		my $find_string = $termnial;
		my $replace_with = join ',' , @split_memers;




		# $third_backbone_copy_newick

		if($count_new_members >= 2)
			{

			if($newick_reformat == 1)# 1 = backbone ID are binomial, convert to genus name only
				{

				if($third_backbone_copy_newick =~ /(.{10}$find_string.{20})/)
					{
					my $context = $1;$context_print++;
				#	if($context_print <= 6){print BCN_LOG "\tfound location in newick to replace:$context\n"};
					};

				# if($count_new_members <= 9){$constreaint = ""}; # !!!!

				$forth_backbone_copy_newick =~ s/([\(\,])$find_string[_][^\:\(\)\,]+/ $1 . "(" . $replace_with . ")" . $constreaint /e;

					# this is the important string replace:
				if($process_backbone_tree_terminal_IDs == 0)
					{
					if($third_backbone_copy_newick =~ s/([\(\,])$find_string(\:)/ $1 . "(" . $replace_with . ")" . $2/e)
						{
						$sucessful_data_insertions++;
						$list_all_inserted_things .= ",$replace_with,";
						}else{$unsucessful_data_insertions++};			
					
					}else{
					if($third_backbone_copy_newick =~ s/([\(\,])$find_string[_][^\:\(\)\,]+/ $1 . "(" . $replace_with . ")" /e)
						{$sucessful_data_insertions++}else{$unsucessful_data_insertions++};			
					};
				}else{
				die "\nnot implemented\n";
				};			


			}else{
			if($find_string eq $replace_with)
				{
				$list_all_inserted_things .= ",$replace_with,";

				}else{

				if($newick_reformat == 1)# 1 = backbone ID are binomial, convert to genus name only
					{
					# remember, higher taxa might have been assigned.
					# and if only one genus is in the new data,
					# it may not neccessary be the same genus as the one in the backbone.
					# hence these might be different here, still need to replcae
					if($process_backbone_tree_terminal_IDs == 0)
						{
						$third_backbone_copy_newick =~ s/([\(\,])$find_string(\:)/ $1 . $replace_with . $2/e;			
						$list_all_inserted_things .= ",$replace_with,";

						}else{
						$third_backbone_copy_newick =~ s/([\(\,])$find_string[_][^\:\(\)\,]+/ $1 . $replace_with /e;			
						};
					}else{
					die "not impleented\n";
					};
				};

			};

		foreach my $mem(@split_memers)
			{
			if($new_members_assigned{$mem} =~ /./){die "\nerror, $mem assigned twice\n"};
			$new_members_assigned{$mem} = $constreaint;
		#	print "$mem contrained to $constreaint\n";
			};

		$doesnt_need_pruning++;
		}else{

		# this array not used in analysis, just list printed to file.
#		push @needs_pruning , $constreaint;
		push @needs_pruning , $termnial; # bugfix 201908
		print GRAFT_LOG "no barcode data for current,CONSTR number:$cosntr_number, termnial:$termnial, taxon_assigned:$constreaint\n";
		print BCN_LOG "terminal $cosntr_number, $termnial, taxon_assigned:$constreaint ... this terminal of backbone will be REMOVED because there was no significant COI data found for it.\n";
	#	print "CONSTR:$constreaint, NO new membres\n";
		};


	########################################################################################################################

	# 20220728: add function to graft source trees to backbone according to overlap of higher taxa.

	if($backbone_terminal_rename2{$termnial} =~ /\w/)
		{
		my $lin5 = $backbone_terminal_rename2{$termnial};

		 print GRAFT_LOG "\tretreived lineage for backbone terminal ($termnial):$lin5\n";
		# retreived lineage for backbone terminal (Trichocera_saltator):Trichocera_saltator;Trichocera;Trichoceridae;Trichoceroidea;
		# this needs to work:
		# retreived lineage for backbone terminal (Tanytarsus_sp):Tanytarsus;

		my $current_details = ""; 
		my $next_highest_ranktax; my $highest_ranktax;my $backbone_rank_skip;
		if($lin5 =~ /(\w+)\;(\w+)\;$/)
			{
			$next_highest_ranktax =  $1; $highest_ranktax =  $2;
			$backbone_rank_skip =0;
			#                              key=taxon             entry= nodeID;count
			# $source_tree_nodes_for_taxon{$taxon_for_node} .= "$current_node;$count_termnials";
			if($source_tree_nodes_for_taxon{$highest_ranktax} =~ /\w/)
				{
				$current_details = $source_tree_nodes_for_taxon{$highest_ranktax};
				}elsif($source_tree_nodes_for_taxon{$next_highest_ranktax} =~ /\w/)
				{
				$current_details = $source_tree_nodes_for_taxon{$next_highest_ranktax};$backbone_rank_skip=1;
				};
			}elsif($lin5 =~ /([A-z][a-z]+)\;$/) # for _sp
			{
			$highest_ranktax =  $1;	
			$current_details = $source_tree_nodes_for_taxon{$highest_ranktax};
			$backbone_rank_skip =0;
			};


			if($current_details =~ /./)
				{

				print GRAFT_LOG " Highest tax ($highest_ranktax) assigneable to terminal ($termnial); details from source trees:$current_details\n";				
				#############################################################
				$current_details =~ s/\t$//;
				my @array_source_nodes_current_tax = split /\t+/, $current_details;
				my $largest_clade_count = 0; my $largest_clade_ID = "";
				foreach my $source_node(@array_source_nodes_current_tax)
					{
					if($source_node =~ /(.+)\;(\d+)/)
						{
						my $nodeID9 = $1; my $count9 = $2; 
						if($count9 > $largest_clade_count)
							{
							$largest_clade_count = $count9;$largest_clade_ID = $nodeID9;
							}
						}else{	
						die "\ncant extract details from string:$source_node\n";
						};
					};

				if($largest_clade_count >= $grafting_cutoff)
					{
					print "\t\tbackbone terminal ($termnial, $constreaint). source node $largest_clade_ID is preferred, it has $largest_clade_count terminals\n";

					if( $sourcetree_subtrees{$largest_clade_ID} =~ /\w/)
						{
						$sourcetree_subtree_retrieved_for_nodeID++;
						my $subtree9 = $sourcetree_subtrees{$largest_clade_ID};
						print GRAFT_LOG "\tbackbone terminal $termnial, assigned higher tax $constreaint / $highest_ranktax / $next_highest_ranktax, rank skip:$backbone_rank_skip, ",
							"will be replaced with source subtree $largest_clade_ID\n\t$subtree9\n";
						if($fourth_backbone_copy_newick =~ s/$termnial/$subtree9/){}else{print "\nwarning 4615.\n"};
						my $cladelabel = $largest_clade_ID; $cladelabel =~ s/INTERNAL_NODE_([\w\d]+).nwk\d+/$1/;
						$backbone_copy_newick5 =~ s/$termnial/$subtree9$constreaint.$cladelabel/;
						$backbone_copy_newick6 =~ s/$termnial/$subtree9$constreaint.$cladelabel/;
						$grafts_made++; $grafted_count += $largest_clade_count;

						# $graft_node_labels{$nodeID}

						}else{
						$sourcetree_subtree_not_retrieved_for_nodeID++;
						};
					};
				#############################################################
				}else{; #if($source_tree_nodes_for_taxon{$highest_ranktax} =~ /\w/)
				print GRAFT_LOG " nothing returned from source tree has for this taxon\n"; 
				};

		#	}else{# if($lin5 =~ /\;(\w+)\;$/)
		#	print GRAFT_LOG "\tnot much lineage info returned.\n";
		#	};

		}else{
		print "NO lineage for backbone terminal ($termnial)\n";
		};



	########################################################################################################################

	};

print "
sourcetree_subtree_retrieved_for_nodeID:$sourcetree_subtree_retrieved_for_nodeID
sourcetree_subtree_not retrieved_for_nodeID:$sourcetree_subtree_not_retrieved_for_nodeID
";
 # taxon replacments with barcode IDs have been made in newick string: $third_backbone_copy_newick

my @new_members_assigned_to_backbine_constraint = keys %new_members_assigned;

my $printthisstring = "
done.
 sucessful_data_insertions:$sucessful_data_insertions++
 unsucessful_data_insertions:$unsucessful_data_insertions++
 new_members_assigned_to_backbone_constraint:$#new_members_assigned_to_backbine_constraint
";

# count doesnt_need_pruning:$doesnt_need_pruning. count need pruning:$#needs_pruning, list:\n@needs_pruning
# these are listed in file XXXXX.needs_pruning

print $printthisstring;
print BCN_LOG $printthisstring;


# close BCN_LOG;
close GRAFT_LOG;


if($graft_source_trees_list =~ /\w/)
	{

open(OUT621, ">backbone_with_grafted_sourcetrees") || die "\nerror 4661\n";
print OUT621 "$fourth_backbone_copy_newick\n";
close OUT621;
open(OUT622, ">backbone_with_grafted_sourcetrees.internal_labelled") || die "\nerror 4661\n";
print OUT622 "$backbone_copy_newick5\n";
close OUT622;
open(OUT6222, ">backbone_with_grafted_sourcetrees.internal_labelled.2nd") || die "\nerror 46611\n";
print OUT6222 "$backbone_copy_newick6\n";
close OUT6222;

	};



 open(TEST_FILE , ">$outfile_prefix.needs_pruning") || die "\nerror 759\n";
foreach my $taxas(@needs_pruning)
	{
	print TEST_FILE ">$taxas\nACTG\n"
	}
close TEST_FILE;



open (OUT69999999 , ">$outfile_prefix.less_basic_constraint_tree") || die "\n\n";
print OUT69999999 "$third_backbone_copy_newick\n";
close OUT69999999;

open (OUT699999998 , ">$outfile_prefix.less_basic_constraint_tree.internal_tax_labeled") || die "\n\n";
print OUT699999998 "$forth_backbone_copy_newick\n";
close OUT699999998;

print "
1 output now printed to file:$outfile_prefix.less_basic_constraint_tree
  also informative is $outfile_prefix.less_basic_constraint_tree.internal_tax_labeled
";


$list_all_inserted_things =~ s/^\,+//;$list_all_inserted_things =~ s/\,+$//;
my @splitit = split /\,+/ , $list_all_inserted_things;
open(FLE, ">$outfile_prefix.list_constrained_members") || die "";
print FLE "@new_members_assigned_to_backbine_constraint\n";
foreach my $thing(@new_members_assigned_to_backbine_constraint) # @splitit)
	{
	print FLE "$thing\t$new_members_assigned{$thing}\n";	
	};
close FLE;


print TESTING "
although script has lots of code from various method implementations, currently only doing one thing.
	check taxa assigned to backbone tree in file:
	$outfile_prefix.less_basic_constraint_tree.backbone_taxa_assigned

	and constraint tree is
	$outfile_prefix.less_basic_constraint_tree
	which will probably require some pruning before doing anything else.
	
	backbone_constraints_newick_LOG list each constraint, and taxa assigned to it


";


};	# sub backbone_constraints2






#####################################################################################################
#
#
#
#####################################################################################################






#
# bottom up. to list individual constraints (ie fasttree format), pseudocode:
#
# recurse through all nodes in root to tip direction, 
# for each node:
#	for each child branch of node:
#		retrive list of taxa assigned to each terminal, for each taxa:
#			retreive list of fasta members assigned to it.
#		set as one grouping of the constraint.
#	set node constraint of two or more groups.




sub print_fasttree_format_relational_constraints
{

my $next = $_[0];
my $sum_branchlength = $_[1]; # from root node
my $node_plot_y = $_[2];


$ft_relational_constraints_subcalls++;
if($ft_relational_constraints_subcalls <= 6){print "ft_relational_constraints_subcalls:$ft_relational_constraints_subcalls\n"};
if($ft_relational_constraints_subcalls =~ /000$/){print "ft_relational_constraints_subcalls:$ft_relational_constraints_subcalls\n"};
# print "\nNEW NODE ($next)\n";


my @non_terminal_child_nodes = ();
my @child_nodes_test = ();
my $all_terminals_from_both = "";
my @taxon_for_child_nodes = ();

my $lowest_taxname_for_node = "";
my $lineage_for_node = "";
$taxa_for_terminals_of_node = "";


if($next =~ /INTERNAL_NODE_\d/) 	# && $all_child_nodes_are_internal == 1)
	{
	}else{#if($next =~ /INTERNAL_NODE_\d/)
	# for tips, you still need to know the lineage. but it wont be a shared lineage.
	# in this case you dont have a list of terminals, theres just one, the curent node
#	$taxa_for_terminals_of_node = get_all_taxa($next);# bug fix 22 oct 2015
	};

if($node_label{$parent_node} =~ /\w/)
	{
	$parent_taxon = $node_label{$parent_node};
	#print "node($next), tax($lowest_taxname_for_node) parent($parent_node) Tax($parent_taxon)\n";
	}else{#if($node_label{$parent_node} =~ /\w/)
	
	}



my $child_nodes = $child_counts{$next};
my @next1 = ();
my @current_node_relational_constraint = ();

# record the taxonomic name(s) given to this node
$nodes{$next}{node_taxonomic_label}=$label_string_to_assign;


my $terminal_members_found = 0;
my @new_collapsed_node;
my $visualization_branch_color;
my %branch_colors;my %branch_color_counts;my %count_barcodes_for_terminals;my %taxname_assigned;

################################################################################################################
for $i( 0 .. $child_nodes )
	{
	# specify new array of the child nodes as normal for recursing
	my $push_node = $nodes{$next}{$i};
	push @next1 , $push_node;	
	$all_terminals_from_this_node =""; # NOTE, this variable is new barcode members, not terminals of original backbone tree.
	my $count_barcodes_this_backbone_child;
	my $constreaint2;

	if( $push_node =~ /INTERNAL_NODE_\d/ ) 
		{
		#############################################
		get_terminals_from_this_node6($push_node,0);#
		#############################################

		# store child internal nodes **NODE_ID** only if it has decendents with barcodes.
		if($all_terminals_from_this_node =~ /[\w\d]/){push @new_collapsed_node , $push_node};

		}else{

		my $constreaint = $store_all_constraints{$push_node};$constreaint2 = $constreaint; # print "the termnial:$push_node the constreaint:$constreaint\n";
		if( $fasta_ids_for_each_taxon{$constreaint} =~ /\S/ )	
			{
			my $new_members_included_in_constraint = $fasta_ids_for_each_taxon{$constreaint};
			$all_terminals_from_this_node .= "	$new_members_included_in_constraint	";

			# store child terminal **Barcode_IDs** if they are present.
			my $barcode_terminals_string = $new_members_included_in_constraint;
			$barcode_terminals_string =~ s/^\t+//;$barcode_terminals_string =~ s/\t+$//;
			my @splitmembers  = split /\t+/, $barcode_terminals_string;
			$count_barcodes_this_backbone_child = scalar @splitmembers;

			if($barcode_terminals_string =~ s/\t+/,/g){$barcode_terminals_string = "(" . $barcode_terminals_string . ")"}; # if there is a comb, these need encasing in parentheses
			push @new_collapsed_node , $barcode_terminals_string;
			foreach my $mem4(@splitmembers){$another_list_of_RelConstrained_barcodeIDs{$mem4} = 1};
			};


		};

	if($all_terminals_from_this_node =~ /[\w\d]/)
		{
		$terminal_members_found++; 
		         $branch_color_counts{$i}++;
		$count_barcodes_for_terminals{$i} = $count_barcodes_this_backbone_child;
		            $taxname_assigned{$i} = $constreaint2;
		push @current_node_relational_constraint , $all_terminals_from_this_node;
		};

	};# for $i( 0 .. $child_nodes )
################################################################################################################




my $find_string = $next;

if($terminal_members_found >= 2) # test is current backbone node has a bipartition in which each contains new barcode members
	{
#	$visualization_branch_color = "black";
#	print "\nFT relational contraint, bipart:\n";
	for my $j(0 .. $#current_node_relational_constraint) # this is list of child internal nodes with backbone decendents
		{
		my $contr98 = $current_node_relational_constraint[$j]; # print "j:$j contr98:$contr98\n";
		$relational_contraints{$next} .= "\t$contr98\tRelationalConstraint$j\t"; # this hash is used for fasttree format file
		};

	# here build up newick string; $pruned_constraint_newick is the main outfile of this script
	my $replace_with = join ',' , @new_collapsed_node;
	if($pruned_constraint_newick =~ s/([\(\,])$find_string([\:\(\)\,])/ $1 . "(" . $replace_with . ")" . $2/e)
		{$sucessful_constraint_insertions_A++};
	if($replace_with =~ /Diadegma/){print "$find_string\tDiadegma\n"}; # INTERNAL_NODE_766	Diadegma


	}elsif( $terminal_members_found == 1 ) 
	{
#	$visualization_branch_color = "gray";
#	print "collapse node\n";		
	my $replace_with = join ',' , @new_collapsed_node;
	if($pruned_constraint_newick =~ s/([\(\,])$find_string([\:\(\)\,])/ $1 .  $replace_with  . $2/e)
		{$sucessful_constraint_insertions_B++};
	if($replace_with =~ /Diadegma/){print "$find_string\tDiadegma\n"};

	}else{
#	$visualization_branch_color = "lightgray";
	};




if(exists($child_counts{$next}))
	{

	for my $index(0 .. $#next1)
		{
		my $test = $next1[$index];




###################################################################################################################################################
#
# 2021-06-06: when dealing with such large datasets, can be rather unclear what has actually happened at this step.
#		So, would be useful to have a visualization. 
#		Will try that here, using tree plotting commands from process_newick.pl
#		At most basic, should be 4 values needed for rectangular plotting:
# 		$sum_branchlength, $y1_proportion , $y2_proportion, $assign_new_x
#
#		1) sum_branchlength (X LEFT) is the easy one, it is input to this subroutine
#		2) y1_proportion requires node_plot_y
#			2a) node_plot_y is 1 for the root node, subsequent nodes are the previous node's assign_new_y
# 			2b) assign_new_y requires 	$index, which is index of child node
#							count_terminals_node_lable
#		count_terminals, takes a single value. calculated in sub record_tree2

# $multipleier 			= 6; # [Circumferal range]

 ######################################
 # X
my $branchlength_to_child;
# my $count_nodes_to_tip_from_current = $nodes_to_tip_from_current{$test} + 1;
# my $remaining = 1 - $sum_branchlength;
my $default_bl = 0.001;
my $bl = $default_bl;
if($branchlengths{$next}{$test} =~ /\d/)
	{$bl = $branchlengths{$next}{$test}};
$branchlength_to_child = $bl;
my $assign_new_x = $sum_branchlength+$branchlength_to_child;
 ######################################
		

 ######################################
 # Y
my $y1_proportion = ($node_plot_y / $count_terminals)+1;	# $y1_proportion *= $multipleier;
my $terminals_from_parent = $count_terminals_node_lable{$next};
my $assign_new_y;
my $current_sum = 0; my $distance;
for my $index_again(0 .. $index)
	{
	my $test_again = $next1[$index_again];
	my $terminals_from_child_again = $count_terminals_node_lable{$test_again};
	unless($terminals_from_child_again =~ /\d/){$terminals_from_child_again = 1};
	my $halfway = $terminals_from_child_again / 2;
	$distance = $current_sum + $halfway;
	#	print "\tindex_again:$index_again of $#next2, count terminals:$terminals_from_child_again\n";
	$current_sum += $terminals_from_child_again;
	};
$delta_y = ($terminals_from_parent / 2) - $distance;
$assign_new_y = $node_plot_y + $delta_y;
my $y2_proportion = ($assign_new_y / $count_terminals)+1;	# $y2_proportion *= $multipleier;
 ######################################

# print "plotting variables 1:$sum_branchlength 2:$assign_new_x 3:$y1_proportion ($node_plot_y / $count_terminals) 4:$y2_proportion ($assign_new_y / $count_terminals)\n";

my $count_new_barcodes;my $taxasign;
if($branch_color_counts{$index} >=1)
	{
	$visualization_branch_color = "black";
	$count_new_barcodes = $count_barcodes_for_terminals{$index};
	$taxasign = $taxname_assigned{$index};
	}else{
	$visualization_branch_color = "gray"
	};

if($count_new_barcodes >= 1)
	{
	$corresponding_barcodes{ $y2_proportion } = $count_new_barcodes;
	$tax_of_Y{ $y2_proportion } = $taxasign;
	};

$rect_branch_color_default = $visualization_branch_color;
$rectangular_plot_branch_width = 2;

# adjust so x is zero to one
my $x_adjust1 = $sum_branchlength / $maximum_x;my $x_adjust2 = $assign_new_x / $maximum_x; 
my $x_adjust3 = $x_adjust2 + 0.15; # terminal labels a bit further out than end of branch

if($y2_proportion >= $plot_y_max){$plot_y_max = $y2_proportion};
if($y1_proportion <= $plot_y_min){$plot_y_min = $y1_proportion};

# parent to one child, vertical:
my $R_command =  "segments(" . "$x_adjust1, $y1_proportion , $x_adjust1, $y2_proportion , col = \"$rect_branch_color_default\", lwd = $rectangular_plot_branch_width)\n";
$draw_tree_R_commands_RECTANGL .= $R_command;
# parent to one child, horizontal
my $R_command =  "segments(" . "$x_adjust1, $y2_proportion , $x_adjust2, $y2_proportion , col = \"$rect_branch_color_default\", lwd = $rectangular_plot_branch_width)\n";
$draw_tree_R_commands_RECTANGL .= $R_command;

if($test =~ /INTERNAL_NODE_\d/) 	# && $all_child_nodes_are_internal == 1)
	{
	}else{
	my $R_command = "text($x_adjust3, $y2_proportion,labels=\"$test\",cex = 0.4, col = \"$rect_branch_color_default\")\n";
	if($terminal_labels_on_R_visualization == 1){$draw_tree_R_commands_RECTANGL .= $R_command};
	};

###################################################################################################################################################









		#####################################
		print_fasttree_format_relational_constraints($test , $sum_branchlength+$branchlength_to_child , $assign_new_y);#	# recurse
		#####################################

		}# 	for my $index(0 .. $#next1)



	};

return();

}


#####################################################################################################
#
#
#
#####################################################################################################



sub get_terminals_from_this_node6 	# 	IS USED (within print_fasttree_format_relational_constraints)
{
my $next = shift;

my $child_nodes = $child_counts{$next};
my @next1 = ();
for $i(0 .. $child_nodes)
	{
	push @next1 , $nodes{$next}{$i};	
	#print "line1423node:$next i:$i child:$nodes{$next}{$i}\n";
	}



for my $index(0 .. $#next1)
	{
	my $test = @next1[$index];#print "test 1430:$test\n";

	if($test =~ /^INTERNAL_NODE_/)
		{

		#####################################
		get_terminals_from_this_node6($test );#	# recurse
		#####################################

		}else{
		$count_terminal_nodes++;

		if($test eq "")
			{
			print "\nnext:($next)\n";
			die "\nerror 4435. looks like sub get_terminals_from_this_node6 was called with a terminal node.\n"
			}



		my $constreaint = $store_all_constraints{$test};$cosntr_number++;
	#	print "termnial:$termnial constreaint:$constreaint\n";

		if($fasta_ids_for_each_taxon{$constreaint} =~ /\S/)	
			{
			my $new_members_included_in_constraint = $fasta_ids_for_each_taxon{$constreaint};
		#	$new_members_included_in_constraint =~ s/^\t+//;
		#	$new_members_included_in_constraint =~ s/\t+$//;
		#	$new_members_included_in_constraint =~ s/\t\t+/	/g;
			$all_terminals_from_this_node .= "	$new_members_included_in_constraint	";

		#	my @split_memers = split /\t+/ , $new_members_included_in_constraint;
		#	my $count_new_members = scalar @split_memers;
			}

			}
	}

return();

} # sub get_terminals_from_this_node6



#####################################################################################################
#
#
#
#####################################################################################################


sub store_tax_heirarchy
{
print "
reading users taxonomic file $taxon_table
";

open(TAXTABLE, $taxon_table) || die "\ncant find taxon table ($taxon_table).\n";
my $tax_table_line_count=0;
while (my $line = <TAXTABLE>)
	{

# 	"$child\t" . 		# child node ID
#	"$current_node\t" . 	# parent node id
#	"$name_string\t" . 	# child tax
#	"$parentname\t" . 	# parent tax
#	"$rank\n"; 		# child rank

	$line =~ s/\n//;$line =~ s/\r//;
	if($line =~ /^(.+)\t(.+)\t(.+)\t(.+)\t(.+)/)
		{
		my $child_ID = $1;my $parent_ID = $2;my $child_tax = $3;my $parent_tax = $4;my $child_rank = $5;
		$rank_hash{$child_rank}++;
		if($tax_table_line_count == 0)
			{
			$root_taxon_name = $parent_tax;
			$starting_node = $parent_ID; 			# root node should be first parent node in the file
			$ncbi_nodes{$parent_ID}{name} = $parent_tax; 	# as names are assigned only for child nodes below, 
			}; 						# the root node needs name assigning by parent


#		my $tax_id = $1;my $parent_tax_id = $2;my $rank = $3;
#		print "tax_id:$tax_id parent_tax_id:$parent_tax_id rank:$rank\n";

		my $current_rankcode = 0;
		foreach my $rank_test("order","suborder","infraorder","series","superfamily","family","subfamily","tribe","subtribe","genus",
			"subgenus","species group","species","subspecies")
				{$current_rankcode++;
				if($child_rank eq $rank_test){$ncbi_nodes{$child_ID}{rank_code} = $current_rankcode}
				};

			$ncbi_nodes{$child_ID}{rank} = $child_rank;
			$ncbi_nodes{$child_ID}{parent} = $parent_ID;
			$ncbi_nodes{$parent_ID}{child_nodes} .= "\t" . $child_ID;
	

	
		$child_tax =~ s/[\(\)\,\[\]\'\#\&\/\:\.\-]/ /g;
		$child_tax =~ s/\s\s+/ /g;$child_tax =~ s/\s+$//;$child_tax =~ s/^\s+//;

		$child_tax =~ s/^Candidatus\s(\w+)$/$1/;

		$ncbi_nodes{$child_ID}{name} = $child_tax;


		$tax_table_line_count++;
		}else{
		print "cant parse line:$line\n";
		};

	};

close TAXTABLE;


print "taxon table has been read. 
	root taxon name:$root_taxon_name
 	node stored:$tax_table_line_count

";
}; # tore_tax_heirarchy


#####################################################################################################
#
#
#
#####################################################################################################



sub build_pruned_constraint_string
{


}


#####################################################################################################
#
#
#
#####################################################################################################







#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


# Called thusly:
# 	$traversal_start_node = $root_node; 
# 	assign_names_to_internal_nodes($traversal_start_node , "New_Root");#

sub assign_names_to_internal_nodes
{
my $current_node = $_[0];my $from_parent = $_[1];
my $count_connections = $child_counts{$current_node};
$count_sub_calls++;

# print "\nsub assign_names\n\tcurrent_node:$current_node from_parent:$from_parent\n";



# defined child nodes, these are all connecting nodes except that from direction of origin in traversal.
my @next1 = ();
my %check_duplication = ();
for my $all_connections(0 .. $count_connections)
	{
	my $connecting_node = $nodes{$current_node}{$all_connections};
	if(exists($check_duplication{$connecting_node})){die "\nduplcaution error:$connecting_node\n"};
	$check_duplication{$connecting_node} = 1;

	if($from_parent eq "New_Root")
		{
		# first node encountered, need all connecting nodes
		push @next1, $connecting_node;
		}else{
		unless($connecting_node eq $from_parent)
			{push @next1, $connecting_node};	
		};

	};

if($current_node =~ /INTERNAL_NODE/)
	{
	$terminals_belonging_to_current_node = "";
	$terminals_belonging_to_current_node_which_have_tax_info = "";
	$max_steps_to_tip = 0;

	###############################################################
	get_terminals_from_this_node_FORCOUNT($current_node , $from_parent , 0);#
	###############################################################

	$terminals_belonging_to_current_node =~ s/(\t)\t+/$1/;
	$terminals_belonging_to_current_node =~ s/^\t+//;$terminals_belonging_to_current_node =~ s/\t+$//;
	my @count_terms_array = split /\t/ , $terminals_belonging_to_current_node;
	my $count_termnials = scalar @count_terms_array; # print "count:$count_termnials ";
	$count_terminals_node_lable{$current_node} = $count_termnials;
	$nodes_to_tip_from_current{$current_node} = $max_steps_to_tip;

	}else{
	$nodes_to_tip_from_current{$current_node} = 0;
	}; # if($current_node =~ /INTERNAL_NODE/)




for my $index(0 .. $#next1)
	{
	my $test = $next1[$index];

	if($test =~ /^INTERNAL_NODE_/)
		{
		#####################################
		assign_names_to_internal_nodes($test , $current_node );#	# recurse
		#####################################
		}
	}

return();

}

#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################





sub get_terminals_from_this_node_FORCOUNT
{
my $current_node = $_[0];my $from_parent = $_[1]; my $steps_to_tip = $_[2];$steps_to_tip++;
my $count_connections = $child_counts{$current_node};

my @next1 = ();
my %check_duplication = ();
for my $all_connections(0 .. $count_connections)
	{
	my $connecting_node = $nodes{$current_node}{$all_connections};
	if(exists($check_duplication{$connecting_node})){die "\nduplcaution error:$connecting_node\n"};
	$check_duplication{$connecting_node} = 1;

	if($from_parent eq "New_Root")
		{
		# first node encountered, need all connecting nodes
		push @next1, $connecting_node;
		}else{
		unless($connecting_node eq $from_parent)
			{push @next1, $connecting_node};	
		};

	};


for my $index(0 .. $#next1)
	{
	my $test = $next1[$index];

	if($test =~ /^INTERNAL_NODE_/)
		{
	#	unless(exists($collapse_nodes{$test}))	
	#		{
		#####################################
		get_terminals_from_this_node_FORCOUNT($test , $current_node , $steps_to_tip);#	# recurse
		#####################################
	#		};
		}elsif($test =~ /[\w\d]/)# 20180307: some suspected errors here, thus require a name
		{

		if($steps_to_tip > $max_steps_to_tip){$max_steps_to_tip = $steps_to_tip};

	#	print "CN:$current_node FP:$from_parent steps_to_tip:$steps_to_tip max:$max_steps_to_tip index:$index test:$test\n";

		# store all terminals:
		$terminals_belonging_to_current_node .= "$test\t";

		# store only terminal with tax info:
		my $test_species = $test;$test_species =~ s/^([A-Z][a-z]+)_.+/$1/;
		my $test_taxonomy = $complete_lineage_for_this_species{$test_species};
		if($test_taxonomy =~ /\w+/)
			{$terminals_belonging_to_current_node_which_have_tax_info .= "$test\t"};

		};
	}




return();

};#sub get_terminals_from_this_node_FORCOUNT


#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



sub get_sum_branchlength
{
my $next = $_[0];
my $sum_branchlength = $_[1]; # from root node
my $node_plot_y = $_[2];

my $child_nodes = $child_counts{$next};
my @next1 = ();
for $i( 0 .. $child_nodes )
	{
	my $push_node = $nodes{$next}{$i};
	push @next1 , $push_node;	
	};

if(exists($child_counts{$next}))
	{
	my $branchlength_to_child;
	my $default_bl = 0.001;
	my $bl = $default_bl;
	if($branchlengths{$next}{$test} =~ /\d/)
		{$bl = $branchlengths{$next}{$test}};
	$branchlength_to_child = $bl;
	my $assign_new_x = $sum_branchlength+$branchlength_to_child;

	if($assign_new_x >= $maximum_x){$maximum_x = $assign_new_x};

	for my $index(0 .. $#next1)
		{
		my $test = $next1[$index];
		get_sum_branchlength($test , $assign_new_x , $assign_new_y);#	# recurse
		}# 	for my $index(0 .. $#next1)
	};

return();



}


#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


sub read_source_trees
{
my @source_tree_array = split /\s+/ , $graft_source_trees_list;

print "\nuser opted to graft source trees to backbone ...\n\treading list of ($#source_tree_array) source trees.\n";

if($graft_source_trees_list =~ /\w/)
	{
	open(GRAFT_LOG, ">grafting_LOG.txt") || die "";
	};

foreach my $file9(@source_tree_array)
	{
	my $source_tree_read = 0;
	$current_sourcetree = $file9;

	open(SOURCETREE, $file9) || die "\nerror cant open source tree file:$file9\n";
	print GRAFT_LOG "\nreading source tree $file9\n";

	while(my $line = <SOURCETREE>)
		{
		$line =~ s/\n//;$line =~ s/\r//;
		if($line =~ /\(/)
			{
			my $current_source_tree = $line;$source_tree_read++;print "\nNEW source tree from file $file9 ... \n";

			#######################
			store_sourcetree($file9, $current_source_tree);#
			#######################

			}
		};
	close SOURCETREE;
	unless($source_tree_read == 1){die "\nerror unexpected tree format:$file9\n"};
	};
print "\nsuccessfully read all source trees.\n";

my @sourcetrees_alltaxa = keys %all_taxa_of_sourcetrees;

print GRAFT_LOG "
successfully read all source trees
count_all_internal_nodes_for_sourcetrees:$count_all_internal_nodes_for_sourcetrees
count_all_nodes for which taxa_found (internal_nodes_for_sourcetrees):$count_all_taxa_found_for_internal_nodes_for_sourcetrees
count taxa over all sourcetrees=$#sourcetrees_alltaxa

";

close GRAFT_LOG;


};

#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################

sub store_sourcetree
{
my $file10 = $_[0];my $current__sourcetree = $_[1];
my $interal_node2 = 0;

# my %sourcetree_nodes = ();
# my %sourcetree_child_counts = ();
my $taxinfo_found_for_terminals=0;my $taxinfo_not_found_for_terminals=0;

while ($current__sourcetree =~ s/\(([^\(\)]+)\)([0-9\.]*)/INTERNAL_NODE_$file10$interal_node2/) # processed sumtrees output from a load of raxml boots
	{
	my $node = $1;my $boot = $2; my $nodeID = "INTERNAL_NODE_$file10$interal_node2"; # print "nodeID:$nodeID node:$node\n";
	my @child_nodes = split /\,/ , $node;#print "\@child_nodes:@child_nodes\n";
	 $sourcetree_child_counts{$nodeID} = $#child_nodes;
	# print "\t\tsourcetree_child_counts, key $nodeID, entry $#child_nodes\n";

	for $i(0 .. $#child_nodes)
		{
		$child_nodes[$i] =~ s/\:(.+)//;
		$sourcetree_nodes{$nodeID}{$i} 			= $child_nodes[$i];
		$sourcetree_nodes{$child_nodes[$i]}{parent} 	= $nodeID;
		# print "\t\twriting to hash sourcetree_nodes, nodeID:$nodeID i:$i child_nodes[i]:$child_nodes[$i]\n";
		my $current_child = $child_nodes[$i];
		unless($child_nodes[$i] =~ /INTERNAL_NODE_/)
			{
		#	if($terminals{$child_nodes[$i]} == 1)
		#		{
		#		if($process_backbone_tree_terminal_IDs == 2){print "\nwarning, you have set script to remove species names, this assumes no conspecifics, otherwise expect crash here.\n"};
		#		die "\nfatal error, found duplicate tip in backbone tree:$child_nodes[$i]\n"
		#		};

			my $genus_name = $child_nodes[$i];$genus_name =~ s/^([A-Z][a-z]+)_.+/$1/;
			my $complete_lineage = $complete_lineage_for_this_species{$child_nodes[$i]};
			unless($complete_lineage =~ /\w/){$complete_lineage= $complete_lineage_for_this_species{$genus_name}};
			unless($child_nodes[$i] =~ /^[\w\_]+$/)
				{
				print "\nerror 100 ($file10). no terminal label:$child_nodes[$i]\n";
				}
			if( $complete_lineage =~ /\w\s/)
				{
				while($complete_lineage =~ s/^([^:]+):(\w+)//)
					{

					my $current_rank = $1;my $current_taxname = $2;
				#	$count_taxa_represented_in_terminals{$current_taxname}++;

					}
				$taxinfo_found_for_terminals++;
				}else{
				$taxinfo_not_found_for_terminals++;

				};

			}#unless($child_nodes[$i] =~ /INTERNAL_NODE_/)
		}
	$sourcetree_root_node = $nodeID; # this is a global and will keep the value in the final loop, so the root can be identified later.
	$interal_node2++;
	}; # while ($newick_string =~ s


print "parsed source tree $file10, it has $interal_node2 nodes
";

print GRAFT_LOG " taxinfo_found_for_terminals:$taxinfo_found_for_terminals
 taxinfo_not_found_for_terminals:$taxinfo_not_found_for_terminals
";

#####################################################################################

$traversal_start_node = $sourcetree_root_node; 

print "\ttraversing current source tree from $traversal_start_node\n";
$sourcetree_nodes_traversed=0;
$taxon_returned_for_node = 0;
$taxon_not_returned_for_node =0;
$string_all_taxa_for_sourcetree = "";
$string_all_taxa_for_sourcetree2 = "";
$complete_source_tree_newick_string  = "($traversal_start_node)";

# called from root here
##########################################################
navigate_sourcetree2($traversal_start_node , "New_Root");
##########################################################

print GRAFT_LOG " print finished traversing, has $sourcetree_nodes_traversed nodes
 taxon_returned_for_node:$taxon_returned_for_node
 taxon_not_returned_for_node:$taxon_not_returned_for_node
 string_all_taxa_for_sourcetree:$string_all_taxa_for_sourcetree
";

open(OUT91 , ">$file10.inttaxlabeled") || die "";
print OUT91 "$complete_source_tree_newick_string\n";
close OUT91;

#######################################################################################################################################
sub navigate_sourcetree2
{
my $current_node = $_[0];my $from_parent = $_[1];$sourcetree_nodes_traversed++;
my $count_connections = $sourcetree_child_counts{$current_node};
# print "\t\tsub navigate_sourcetree2  current_node:$current_node count_connections:$count_connections\n";


my @next1 = ();
my %check_duplication = ();
for my $all_connections(0 .. $count_connections)
	{
	my $connecting_node = $sourcetree_nodes{$current_node}{$all_connections};
#	print "\t\t\tconnecting_node:$connecting_node\n";

#	if(exists($check_duplication{$connecting_node})){die "\nduplcaution error:$connecting_node\n"};
#	$check_duplication{$connecting_node} = 1;

	if($from_parent eq "New_Root")
		{
		# first node encountered, need all connecting nodes
		push @next1, $connecting_node; # print "\t\t\tall nodes of root\n";
		}else{
		unless($connecting_node eq $from_parent)
			{
			push @next1, $connecting_node;
			# print "\t\t\talso connecting node $connecting_node\n";
			};	
		};

	};

$taxon_for_node = "";
if($current_node =~ /INTERNAL_NODE/)
	{

	$terminals_belonging_to_current_node = "";
	$source_tree_newick_string = "($current_node)";

	##########################################
	get_terminals_from_node($current_node , $from_parent);#
	##########################################

	# 20220823: seem to have an unneccessary pair of parentheses later:
	$source_tree_newick_string =~ s/^\(//;$source_tree_newick_string =~ s/\)$//;
	$sourcetree_subtrees{$current_node} = $source_tree_newick_string;

	$terminals_belonging_to_current_node =~ s/(\t)\t+/$1/;
	$terminals_belonging_to_current_node =~ s/^\t+//;$terminals_belonging_to_current_node =~ s/\t+$//;
	my @count_terms_array = split /\t/ , $terminals_belonging_to_current_node;
	my $count_termnials = scalar @count_terms_array; 
#	print "\nretrieved all terminals from node $current_node, count:$count_termnials, ($terminals_belonging_to_current_node)\n";

				####################################################
	$taxon_for_node = 	get_shared_taxa($terminals_belonging_to_current_node);#
				####################################################

#	print "\tshared taxon:$taxon_for_node\n";
	# TO DO: need to look into printing which terminals lineages cannot be retreived for

	$count_all_internal_nodes_for_sourcetrees++;
	if($taxon_for_node =~ /\w\w\w/)
		{
		$taxon_returned_for_node++;
		$count_all_taxa_found_for_internal_nodes_for_sourcetrees++;

		# for some tax names there will be lots of nodes assigned it in a given source tree, 
		# dont add to list if already have a larger one.
		# i thinl there is another checl for size later, across all source trees 

	#	print "\nchecking whether to add this to list, current_sourcetree:$current_sourcetree taxon_for_node:$taxon_for_node count_termnials:$count_termnials\n";
		if($count_termnials > $biggest_clade_for_taxon{$current_sourcetree}{$taxon_for_node} && 
			$count_termnials >= $grafting_cutoff ) # dont want to store clades which wont be grafted, or they might interfeer with clades that might be (eg Simuliinae for which there is tiny clade in one tree for higher rank that stops check for subfam)
			{
			$source_tree_nodes_for_taxon{$taxon_for_node} .= "$current_node;$count_termnials\t";
			$biggest_clade_for_taxon {$current_sourcetree}{$taxon_for_node} = $count_termnials;
			print GRAFT_LOG " for source tree $current_sourcetree, storing taxon $taxon_for_node for node $current_node of count $count_termnials decendents\n";
			# print "\tyes\n";
			}else{
			print GRAFT_LOG " already have this taxon recorded for a larger clade $biggest_clade_for_taxon{$current_sourcetree}{$taxon_for_node}\n";
			
			};

		unless($string_all_taxa_for_sourcetree =~ / $taxon_for_node /){$string_all_taxa_for_sourcetree .= " $taxon_for_node "};

		$all_taxa_of_sourcetrees{$taxon_for_node}=1;
		}else{$taxon_not_returned_for_node++};


	}else{
#	$nodes_to_tip_from_current{$current_node} = 0;
	}; # if($current_node =~ /INTERNAL_NODE/)



#################################################################
# would be useful output to have source trees with internal tax labels
my $childnodes55 = join(',' , @next1);
my $replacement_string55 = "(" . $childnodes55 . ")$taxon_for_node";
$complete_source_tree_newick_string =~ s/$current_node/$replacement_string55/;
print GRAFT_LOG "\tassigned taxon $taxon_for_node to node $current_node\n";

#################################################################



for my $index(0 .. $#next1)
	{
	my $test = $next1[$index]; # print "\t\t\trecursing index:$index test:$test\n";

	if($test =~ /^INTERNAL_NODE_/)
		{
		#####################################
		navigate_sourcetree2($test , $current_node );#	# recurse
		#####################################
		}
	}

return();

} # sub navigate_sourcetree2
#######################################################################################################################################







#####################################################################################


}; # sub store_sourcetree


#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################





sub get_terminals_from_node
{
my $current_node = $_[0];my $from_parent = $_[1];
my $count_connections = $sourcetree_child_counts{$current_node};
# print "\t\t\tsub get_terminals_from_node current_node:$current_node count_connections:$count_connections\n";


my @next1 = ();
my %check_duplication = ();
for my $all_connections(0 .. $count_connections)
	{
	my $connecting_node = $sourcetree_nodes{$current_node}{$all_connections};
	if($from_parent eq "New_Root")
		{
		# first node encountered, need all connecting nodes
		push @next1, $connecting_node; # print "\t\t\tall nodes of root\n";
		}else{
		unless($connecting_node eq $from_parent)
			{
			push @next1, $connecting_node;
			# print "\t\t\talso connecting node $connecting_node\n";
			};	
		};

	};

if($current_node =~ /INTERNAL_NODE/)
	{
	}else{
	};





my $childnodes = join(',' , @next1);
my $replacement_string = "(" . $childnodes . ")";
$source_tree_newick_string =~ s/$current_node/$replacement_string/;
# print "replacing $current_node with $replacement_string (child nodes array:@next1)\n";


for my $index(0 .. $#next1)
	{
	my $test = $next1[$index]; # print "\t\t\trecursing index:$index test:$test\n";

	if($test =~ /^INTERNAL_NODE_/)
		{
		#####################################
		get_terminals_from_node($test , $current_node );#	# recurse
		#####################################
		}else{
		$terminals_belonging_to_current_node .= "$test\t";
		
		}
	}

return();

}


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





sub record_tree_graft
{


my $tree1= "";
open(IN, $treefile) || die "\n\nerror 1408 cant open file $treefile\n\n";
while (my $line = <IN>)
	{
	if($line =~ /(.+\(.+)/){$tree1 = $1}
	}
close(IN);

print "\nsub record_tree_graft. looking at tree:$treefile\n";

$tree1 =~ s/ //g;

$second_backbone_copy_newick = $tree1;
$third_backbone_copy_newick = $tree1;  # for inserting barcode IDs
$fourth_backbone_copy_newick = $tree1; # for grafting source trees

# this one will have internal taxonomic node labesl:
$forth_backbone_copy_newick = $tree1;

# thus , rm supprto node labesl:
$forth_backbone_copy_newick =~ s/(\))[0-9\.]+/$1/g;

# grafting with internal labels
$backbone_copy_newick5 = $tree1;
$backbone_copy_newick5 =~ s/(\))[0-9\.]+/$1/g;

# grafting with internal labels both current and prev iteration
$backbone_copy_newick6 = $tree1;
$backbone_copy_newick6 =~ s/(\))[0-9\.]+/$1/g;

# print "\n$tree1\n";die "";

$tree_parse = 1; 	# assuming sensible tree.

if($plot_constraints == 1){$tree_parse = 1}; # user more interested in the R plot which only works with no lengths

if($tree_parse ==1)
	{
	# remove branchlengths, scientific notation, incl negative values for distance trees. example: -8.906e-05
#	$tree1 =~ s/\:\-*\d+\.\d+[eE]\-\d+([\(\)\,])/$1/g; 
	# remove regular branchlengths: 0.02048
	$tree1 =~ s/\:\-*\d+\.\d+//g; 
	# remove 0 length branchlengths
	$tree1 =~ s/\:\d+//g;
	}else{
	print "script set to not remove all formats of branch length, this assumes sensible newick string\n";
	};


my $newick_string 	= $tree1;
my $interal_node	= 0;
$backbone_copy_newick = $tree1;






#####################################################################################

# new newick parser ....

#while ($newick_string =~ s/\(([^\(\)]+)\)(\d*)/INTERNAL_NODE_$interal_node/) # rooted bipart raxml tree, the 2 last nodes will have no boot values
#	{

while ($newick_string =~ s/\(([^\(\)]+)\)([^\)\(\,]*)/INTERNAL_NODE_$interal_node/) # processed sumtrees output from a load of raxml boots
	{

	my $node = $1;my $nodelable = $2; my $nodeID = "INTERNAL_NODE_$interal_node"; # print "nodeID:$nodeID node:$node\n";


		# seems raxml constraint needs to be unrooted
	if($newick_string =~ s/^\((INTERNAL_NODE_\d+)\:[\d\.]+(\,[A-Za-z_]+\:[\d\.]+)\)(\;)/$1$3/)
		{
		$node .= $2;
		print "\nwarning, since raxml constraint needs to be unrooted, " ,
			"tried to unroot your rooted tree. may or may not work.\n\troot node now:$node\n";
		};

	if($nodelable =~ /\w/)
		{
	#	print "\nnodelable:$nodelable\n";
		$graft_node_labels{$nodeID} = $nodelable;
	#	$nodes{$nodeID}{support} = $boot
		}else{
	#	$nodes{$nodeID}{support} = 100
		};#print "\tnodeID:$nodeID boot:$boot\n";

	my @child_nodes = split /\,/ , $node;#print "\@child_nodes:@child_nodes\n";
	$child_counts{$nodeID} = $#child_nodes;

	if($interal_node >= 6)
		{}elsif($interal_node == 5)
		{print "........\n.......\n"}else{
		print "nodeID:$nodeID\n\t$#child_nodes child_nodes:@child_nodes\n";
		};


	for $i(0 .. $#child_nodes)
		{
		#print "$child_nodes[$i]\n";
		my $bl = "NA"; 	
		if($child_nodes[$i] =~ /\:(.+)/)	# branchlength found
			{$bl = $1};
		$child_nodes[$i] =~ s/\:(.+)//;
		# optional processing of tip labels
		unless($child_nodes[$i] =~ /INTERNAL_NODE_/)
			{
			if($remove_accessions_from_reference_species == 1)
				{$child_nodes[$i] =~ s/^([A-Z][a-z]+_[a-z]+)_.+$/$1/}
			if($process_backbone_tree_terminal_IDs == 1)
				{$child_nodes[$i] =~ s/.+_([A-Z][a-z]+_[a-z]+)$/$1/}
			if($process_backbone_tree_terminal_IDs == 2)
				{
				unless($child_nodes[$i] =~ s/^([A-Z][a-z]+)_.+$/$1/){die "\nERROR, process_backbone_tree_terminal_IDs == 2, however cant remove species from terminal:$child_nodes[$i]\n"};
				
				}
			};


		# record each decendent of the current node, in a hash nodes{ the current node }{ child number }
		$nodes{$nodeID}{$i} 			= $child_nodes[$i];
		$nodes{$child_nodes[$i]}{parent} 	= $nodeID;
		$nodes{$child_nodes[$i]}{branchlength} 	= $bl;


		# record length of branch between the current two nodes,
		# store in both directions
		my $current_child = $child_nodes[$i];
		$branchlengths{$nodeID}{$current_child} = $bl; # print "storing bl ($bl) to $nodeID / $current_child\n";
		$branchlengths{$current_child}{$nodeID} = $bl;


		# and record whether the current node is a terminal

		unless($child_nodes[$i] =~ /INTERNAL_NODE_/)
			{
			if($terminals{$child_nodes[$i]} == 1)
				{
				if($process_backbone_tree_terminal_IDs == 2){print "\nwarning, you have set script to remove species names, this assumes no conspecifics, otherwise expect crash here.\n"};
				die "\nfatal error, found duplicate tip in backbone tree:$child_nodes[$i]\n"
				};

			if($child_nodes[$i] =~ /[A-Z][a-z]+_[A-Z][a-z]+_[A-Z][a-z]+_/){$poorly_labelled_terminals++};
			if($poorly_labelled_terminals >= 20)
				{die "\nERROR. your terminal labels (e.g. $child_nodes[$i]) cannot be parsed. should be no more than Genus_species. quitting.\n\n "};
 	
			$terminals{$child_nodes[$i]} =1;$count_the_terminal_nodes++;
			my $genus_name = $child_nodes[$i];$genus_name =~ s/^([A-Z][a-z]+)_.+/$1/;
			my $complete_lineage = $complete_lineage_for_this_species{$child_nodes[$i]};
			unless($complete_lineage =~ /\w/){$complete_lineage= $complete_lineage_for_this_species{$genus_name}};
			unless($child_nodes[$i] =~ /^[\w\_]+$/)
				{
				print "\nerror 100. no terminal label:$child_nodes[$i]\n";
				# die "";
				}
			if( $complete_lineage =~ /\w\s/)
				{
				$complete_lineages_retreived++;
				while($complete_lineage =~ s/^([^:]+):(\w+)//)
					{
					my $current_rank = $1;my $current_taxname = $2;
				#	print "\tcurrent_rank:($current_rank) current_taxname:($current_taxname)\n";
				#	current_rank:( family) current_taxname:(Bothrideridae)
					$count_taxa_represented_in_terminals{$current_taxname}++;
				#	$tax_is_of_this_rank
					}

				}else{
				$errors_printed++;
				if($errors_printed < 6)
				{print "\nERROR (6075).  no complete lineage retrieved for genus $genus_name\n" , 
				"\t$complete_lineages_retreived have sucessfully been retrvied previously\n" , 
				"\ttip: check you are not using outdated NCBI taxonomy database; check taxon specified in taxon table ($root_taxon_name) matches those in your tree\n";
				print "check its in there with command grep \"$genus_name\" names.dmp\n";
				}elsif($errors_printed == 6){print "..... not printing further errors of this type.\n"};

				if($genus_name =~ /^Beetles$|^Flies$|^Moths$|^Wasps$/)
					{
					}else{
				#	die ""
					};
				# die "\tquitting.\n\n";
				};



			}#unless($child_nodes[$i] =~ /INTERNAL_NODE_/)

		}

	#print "node:$interal_node\n\tchild nodes:@child_nodes\n";
	$root_node = $nodeID; # this is a global and will keep the value in the final loop, so the root can be identified later.
	$interal_node++;$print_internal_nodes++;
#	if(length($newick_string)<= 500){print "newick_string:$newick_string\n"};

	}; # while ($newick_string =~ s

#####################################################################################

print "
taxonomy not retrived for $errors_printed terminals
";

unless($complete_lineages_retreived >= 1)
	{
	print "\n\n\nerror, your taxon table does not contain any of the things in your phylogeny. " , 
		"please check root tax number input to taxon_table.pl. quitting.\n\n\n";
		die "";
	};


 print "remaining of newick_string :$newick_string\n";

my @terminal_array = keys %terminals; @terminal_array = sort @terminal_array;
# store single value, used later:
$count_terminals = scalar @terminal_array;



unless($interal_node >= 2){die "\nerror reading your phylogeny.\n"}


}#sub record_tree_graft




#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################





