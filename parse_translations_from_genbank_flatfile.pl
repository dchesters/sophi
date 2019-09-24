#!/usr/bin/perl



# 
# 
# 
# 
# 
# 
# 
# 2015/10/06: inserted code for organizing input seqs to output files by gene name
# 	some entries have duplicates, e.g. COX1 in twice: Watasenia scintillans; NC_007893.1  GI:89337394
# 	onyl print one in these cases.
# 2016/04/09: all settings now given in command, none written into script itself
# 2016/05/20: bugfix, species names were being truncated one character
# 2017/05/18: ncbi discontinued GI numbers, parse accession instead
# 2017/05/19: corrected some problem names
# 2018/05/04: outputs excised dna sequence for each feature
# 
# 
# 
# 
# 
# 
# 
# 
# 
#########################################################################################################






my $arg_string  = join ' ', @ARGV;


# $id_format = 5; # 1=protein id only, 2=gi then protein id
	# 3= new unique number to each seq, 4= species name etc 
	# 5 is standardized protein name 
# $organize_outfiles_by_genename = 1;
# $dont_print_duplicate_seq = 0;
# $prefered_field_for_parsing_gene_labels = 1;# 1= gene < note < product. original setup, used a lot for insect mtgenomes
	# 1 still works better for mitogenomes, just been trying Dikarya, if using option 2 multiple sequences are found for many names, when should be just one each
	# 2 = note < product < gene. seems better while using chloroplasts 









$sequence_count=0;

# problems, protein ID in features are repeated within genomes






# read a curated list of gene names and alternative names, for the purpose of congruence testing later:
%gene_synonyms;

######################
read_gene_synonyms();#
######################




my $date = localtime time;$month = $date;
$month=~ s/\w+\s+(\w+)\s+\d+\s+\S+\s+(\d+)$/$1$2/;
#$month= "Apr2013";


#####################################################################################################################################

$search_for_accessions = 0;
%accessions_to_look_for;

#####################################
read_command_arguments($arg_string);#
#####################################


my %sp_rm = ();
$countnumberprinted=0;
$counttotal=0;
$inv_number;
my $number_of_this_gene 	= 0;
$zipped_infileTEST ="";




############################
parse_genbank_flatfiles2();#
############################



print "\n\nscript completion with no major error. bye.\n\n";
exit;














##############################################################################################




sub parse_genbank_flatfiles2
{

	####################################
	parse_genbank_flatfileNEW($user_infile);#
	####################################


}




###########################################################################

sub parse_genbank_flatfileNEW
{
my $inputfile = shift;
print "\n\nsub parse_genbank_flatfileNEW\n\n";




open(IN2, $inputfile) || die "cant open infile:$inputfile\n";
print "opened $inputfile\n";

my $fileasstring = "";

while(my $line = <IN2>)
	{
	$fileasstring .= $line;
	}
close(IN2);

my @file_as_array = split(/\n\/\/\n/ , $fileasstring);	# problem here, crashed with a big file, 
							# not sure if its just the size of the file or some feature of.
							# perl error message was 'split loop'


open(USEROUT , ">$user_outfile") || die "";
open(FEATURE_DEBUG, ">DNA_feature_parse_log" ) || die "";
print FEATURE_DEBUG "accession\tspeciesID\tproteinID\tproduct\tfeature_start\tfeature_end\tcomplement\textracted_feature\n";


foreach my $entry(@file_as_array)
	{
	#print "\n$inputfile. from all files, $countnumberprinted were printed out of $counttotal\n";
	#if ($entry =~ /Z93710/i){print $entry}

#if($entry =~ /KC880663|JQ024878|HQ266313/){

%store_proteinIDs = ();

	#########################
	parse_this_entry($entry);#	
	#########################

my @keys3 = keys %store_proteinIDs;
# print "$#keys3 proteins in entry\n";

#die};

	}

close USEROUT;
close FEATURE_DEBUG;

}





###########################################################################



sub parse_this_entry
{
my $current_entry = shift;
my $current_entry_copy = $current_entry;
$counttotal++;




#######		ACCESSIONS

my $accession = "UnknownAccesion";
if($current_entry =~ /ACCESSION\s+(\S+)\n/)
	{
	# single values for accession found
	$accession = $1
	}else{
	# multiple values for accession, first one usually matches with that found on description
	#DEFINITION  Caenorhabditis elegans DNA for 32-kDa galectin, complete cds.
	#ACCESSION   AB000802 D85885 D85886 D85887
	#VERSION     AB000802.1  GI:2564035
	if($current_entry =~ /ACCESSION\s+(\S+)\s(.+)\n/)
		{
		$accession = $1;
		#print "$1 $2\n"
		}else{
		if($current_entry =~ /\w/)
			{print "ERROR, no accession\n$current_entry\n"}	
			}
	}

#VERSION     U49845.1  GI:1293613
my $gi_number = "UnknownGI";
if($current_entry =~ /VERSION\s+\S+\s+GI\:(\d+)\n/)
	{
	# single values for accession found
	$gi_number = $1
	}else{
	 "\nerror 177 NO GI\n$current_entry\nNO GI\n"
	}

# print "\n\nGI_number:$gi_number\n";



#######		TAXID

my $current_taxid = "UnknownTaxID";
my $taxid_present = 0;

	if($current_entry =~ /\/db_xref.+taxon.(\d+)\"/)
		{
		$current_taxid = $1;$taxid_present = 1;
	#	}elsif{
	#	$current_entry =~ /db_xref="BOLD:CNCBR694-09.COI-5P"/ # sometime bold id is assigned to ncbi taxon field. can be ignored, no use to anyone
		}else{
		print LOG2 "\n\nERROR, no tax id\n\n$current_entry\n\n";
		if($current_entry =~ /\w/)
			{
			print "\nERROR, no tax id\n$current_entry\n";
			};
		};

#                     VRVWFCNRRQKEKRINPSLDSPTGADDDESSYMMH"
#ORIGIN      chromosome 2L, map position 3Arabidopsis thaliana|Glycine max|Medicago truncatula|Oryza sativa|Populus trichocarpa|Sorghum bicolor|Vitis vinifera|Zea mays3F polytene chromosome band.
#        1 agtttcacac gagacacgaa cgcgtccaca cgaaagcggc gcttccaaga aatcagaata

if($current_entry =~ /[\n\r]ORIGIN.+[\n\r]\s+1/){$current_entry =~ s/([\n\r]ORIGIN).+([\n\r]\s+1)/$1$2/
	}else{
#	die "no origin:$current_entry\n"
	}


my $current_organism = "";
my $current_taxstring = "";
my $genome = 0;
my $bold = 0;




#######		ORGANISM

# problem labels encountered:
#  ORGANISM  Anopheles culicifacies B
#  ORGANISM  Anopheles dirus A
#  ORGANISM  Anopheles farauti No. 4
#  ORGANISM  Bemisia tabaci complex sp. Asia I
#  ORGANISM  Gasteruption sp. M19
#  ORGANISM  Homolobus sp. QL-2013
#  ORGANISM  Ichneutes sp. QL-2013
#  ORGANISM  Idris sp. MM-2013
#  ORGANISM  Megalyra sp. MM-2014
#  ORGANISM  Neotermes sp. A TB-2014
#  ORGANISM  Pambolus sp. QL-2013



	if($current_entry =~/ORGANISM\s*(.+)\n/)
		{
		$current_organism = $1;

		$current_organism =~ s/^\s*([A-Z][a-z]+ [a-z][a-z][a-z]+) [A-Z]\s*$/$1/;#  ORGANISM  Anopheles culicifacies B
		$current_organism =~ s/^\s*([A-Z][a-z]+ [a-z][a-z][a-z]+) No\. \d+\s*$/$1/;#  #  ORGANISM  Anopheles farauti No. 4


		$current_taxstring = $current_organism;
		$current_taxstring =~ s/^\s//g;$current_taxstring =~ s/\s$//g;
		$current_taxstring =~ s/\s/_/g;
#print "current_taxstring:$current_taxstring\n";
	#	if($current_organism=~/Drosophila melanogaster|Bombyx mori|Tribolium castaneum|Apis mellifera|Anopheles gambiae/i)
	#	if($current_organism =~ /Acyrthosiphon pisum|Aedes aegypti|Anopheles gambiae|Apis mellifera|Bombyx mori|Culex quinquefasciatus|Drosophila|Drosophila melanogaster|Drosophila pseudoobscura|Pediculus humanus|Tribolium castaneum/i)
		#if($current_organism =~ /Arabidopsis thaliana|Glycine max|Medicago truncatula|Oryza sativa|Populus trichocarpa|Sorghum bicolor|Vitis vinifera|Zea mays/i || 
		#	$current_organism =~ /Acyrthosiphon pisum|Aedes aegypti|Anopheles gambiae|Apis mellifera|Bombyx mori|Culex quinquefasciatus|Drosophila|Drosophila melanogaster|Drosophila pseudoobscura|Pediculus humanus|Tribolium castaneum/i
		#	)
		if($current_organism=~/Homo sapiens/i)
			{$genome=1;
		#	print "\tgenome\n"
			}

		if($current_organism =~ /BOLD\:/){$bold = 1
			#;$current_organism =~ s/BOLD\:/BOLD/;
			}

		}

	if($current_organism =~/BOLD/ && $bold == 0){die "error BOLD:$current_entry"}




#######		SPECIES ID

my $current_sp_id = "NA";

	if($current_entry =~/\/organism\=\"([^\"]+)\"/)
		{#/organism="Catenula sp."
		$current_sp_id = $1;$current_sp_id =~ s/[\n\r\.]//g;#$current_sp_id =~ s/^([a-zA-Z\s]+)[\:\,].+/$1/g;
		};







#  *************************************************

### PROTEIN ... and more

my $protein_seq = "";
if($current_entry =~ /\w/)
	{		 ###################################################
	$protein_seq = parse_protein_seqs_in_this_entry($current_entry); # parses multiple features, if present
	}		 ###################################################

# variable $protein_seq is currently inactive, though earlier seems was concatenate of all translations

#  *************************************************


						if($protein_seq =~ /\w/ && $parse_protein == 1)
							{
							open(PROTEIN_OUT , ">>prot") || die "\n836\n";
							print PROTEIN_OUT ">$gi_number" , "_$current_taxid" , "\n$protein_seq\n";
							close PROTEIN_OUT;
              						$number_of_entries_with_protein++;
							print "\nprinting protein seq\n\n";
							}else{
							# print "\nnot printing prot seq!\n";
							}
						# $current_taxstring is the organism name:
						# print MAIN_OUT ">$current_taxstring" , "_$accession" ,"\n$current_dna\n";	
						# gradualy switching to GI, easier to parse than accession
						# you prob want this one for invert mtgenomes:
						 print MAIN_OUT ">$current_taxstring" , "_$gi_number" ,"\n$current_dna\n";	


				
			

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
\n
     ***** parse translations from genbank flatfiles *****
		  
\n";



if($arguments =~ /-accessions_list\s+(\S+)/)
	{
	my $accessions_from_this_file = $1;
	$search_for_accessions =1;
	print "looking for certain members in the .gb files, " . 
	"only printing those present in your file ($accessions_from_this_file)\n\n";

	open(IN, $accessions_from_this_file) || die "\nerror 852\n";
	while(my $line = <IN>)
		{
		if($line =~ />[A-Z][a-z]+_[a-z]+_(.+)$/)
			{
			my $ac = $1;$accessions_to_look_for{$ac}=1; #print "accession:$ac\n";
			}else{
			if($line =~ />/)
				{#print "$line"
				}
			};
		};
	close IN;

	my @countthese = keys %accessions_to_look_for;
	print "recorded $#countthese accessions from your file\n";

	}else{
$search_for_accessions =0;
	}

if($arguments =~ /-in\s+(\S+)/)
	{
	$user_infile = $1;
	}else{die "\ncommand syntax error 1\n"}
if($arguments =~ /-out\s+(\S+)/)
	{
	$user_outfile = $1;
	}else{die "\ncommand syntax error 2\n"}

if($arguments =~ /-out_suffx\s+(\S+)/)
	{
	$out_suffx = $1;
	}else{die "\ncommand syntax error 3\n"}


# $id_format = 5; # 1=protein id only, 2=gi then protein id
	# 3= new unique number to each seq, 4= species name etc 
	# 5 is standardized protein name 
if($arguments =~ /-id_format\s+(\d)/)
	{
	$id_format = $1;
	}else{
	$id_format =5;
	print "\nuser did not give format for IDs in fasta outputs (-id_format [1-5]), using default setting 5 (standardized protein name)\n";
	};

# $organize_outfiles_by_genename = 1;
if($arguments =~ /-outfile_by_annotation/)
	{
	$organize_outfiles_by_genename = 1;
	}else{
	$organize_outfiles_by_genename = 0;
	};

# $dont_print_duplicate_seq = 0;
if($arguments =~ /-filter_duplicates/)
	{
	$dont_print_duplicate_seq = 1;
	}else{
	$dont_print_duplicate_seq = 0;
	};

# $prefered_field_for_parsing_gene_labels = 1;# 1= gene < note < product. original setup, used a lot for insect mtgenomes
	# 1 still works better for mitogenomes, just been trying Dikarya, if using option 2 multiple sequences are found for many names, when should be just one each
	# 2 = note < product < gene. seems better while using chloroplasts 
if($arguments =~ /-preferred_annotation_field\s+(\d)/)
	{
	$prefered_field_for_parsing_gene_labels = $1;
	}else{
$prefered_field_for_parsing_gene_labels = 1;
	};



print "

SETTINGS:
 accession list supplied:$search_for_accessions
 infile:$user_infile
 outfile:$user_outfile
 format for IDs:$id_format
 outfile_by_annotation:$organize_outfiles_by_genename
 filter_duplicates:$dont_print_duplicate_seq
 preferred_annotation_field:$prefered_field_for_parsing_gene_labels

";

#print "\n
#";



}#sub read_command_arguments



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################

sub parse_protein_seqs_in_this_entry
{
my $gb_entry = shift; # print "gb_entry:$gb_entry\n"; 

my $gi_number2 = "UnknownGI";


#if($gb_entry =~ /VERSION\s+\S+\s+GI\:(\d+)\n/)
#	{
#	$gi_number2 = $1
#	}else{
#	open(ERROR_LOG , ">ptfgf_ERRORLOG") || die "\n\n";
#	print ERROR_LOG "error 388 NO GI \n$gb_entry \nerror 388 NO GI";
#	close ERROR_LOG;
#	die "\nerror 388 NO GI \n$gb_entry \nerror 388 NO GI, see LOG file ptfgf_ERRORLOG
# prev this was caused by a truncated ncbi entry, open file $user_infile and delete this
#	\n"
#	}

# 2017 may, ncbi discontinued gi numbers ..

# my $accession = "UnknownAccesion";
if($gb_entry =~ /ACCESSION\s+(\S+)\n/)
	{
	# single values for accession found
	$gi_number2 = $1
	}else{
	# multiple values for accession, first one usually matches with that found on description
	#DEFINITION  Caenorhabditis elegans DNA for 32-kDa galectin, complete cds.
	#ACCESSION   AB000802 D85885 D85886 D85887
	#VERSION     AB000802.1  GI:2564035
	if($gb_entry =~ /ACCESSION\s+(\S+)\s(.+)\n/)
		{
		$gi_number2 = $1
		#print "$1 $2\n"
		}else{
		if($gb_entry =~ /\w/)
			{print "ERROR, no accession\n$gb_entry\n"}	
			}
	}




my $current_sp_id = "NA";
if($gb_entry =~/ORGANISM\s*(.+)\n/)
	{
	$current_sp_id = $1;# print "A current_sp_id:$current_sp_id\n";
	$current_sp_id =~ s/^\s+//;
	$current_sp_id =~ s/\s+$//;$current_sp_id =~ s/\s/_/g;
	$current_sp_id =~ s/^([A-Z][a-z]+_[a-z]+)[^a-zA-Z]+/$1/;# bugfix 20160520
	$current_sp_id =~ s/[^\_\|a-z0-9]//ig;
	# print "\tB current_sp_id:$current_sp_id\n";
	}elsif($gb_entry =~/\/organism\=\"([^\"]+)\"/)
	{#/organism="Catenula sp."
	$current_sp_id = $1;$current_sp_id =~ s/[\n\r\.]//g;#$current_sp_id =~ s/^([a-zA-Z\s]+)[\:\,].+/$1/g;
	$current_sp_id =~ s/^([A-Z][a-z]+_[a-z]+)[^a-zA-Z]+/$1/;
	$current_sp_id =~ s/[^\_\|a-z0-9]//ig;
	}else{print "\nwarning cant pares species.\n"};





my @split_gb_entry = split /\s+CDS\s+|\s+tRNA\s+|\s+rRNA\s+|\s+misc_feature\s+/ , $gb_entry;# 



#####         DNA

my $current_dna = "";
$gb_entry =~ s/\012\015?|\015\012?//g;
$gb_entry =~ s/(ORIGIN\s+)\d bp upstream of [a-zA-Z]{1,10} site\./$1/;
if($gb_entry =~/^.+ORIGIN(.+)$/)
	{
	$current_dna = $1; # $current_entryCOPY2 = $current_entry;
	$gb_entry =~ s/^.+ORIGIN//;
	if($gb_entry =~/^.+ORIGIN(.+)$/)
		{print"\n\nERROR, 2 * ORIGIN?\n\n$accession\n$current_entry_copy";
	 	die "\n\nerror 423\n"}
	$current_dna =~ s/\d//g;$current_dna =~ s/[\s ]//g;
	};
# $gb_entry not used again 










my $concatenated_protein_seqeucnes_for_this_entry = "";
my $number_of_features_found=0;
my $features_without_translation=0;
for  my $k(1 .. $#split_gb_entry)
	{
	my $gb_feature = $split_gb_entry[$k];	# print "\n\nFEATURE:($gb_feature)\n";

	my $proteinID = "NA";#                     /protein_id="WP_025368665.1"
	if($gb_feature =~ /\/protein_id\=\"([^\"]+)/)
		{$proteinID = $1};

 #                    /db_xref="GI:644811070"



	$gb_feature =~ s/\s\s+gene\s\s\s+.+$//s;

# presumably above regex was supposed to trim off:
#                     EDDFFSIKLHVKELEDASKDFFSLRLKKSINSFTKK"
#     gene            640877..641212
#                     /locus_tag="BUMPG002_RS03095"

# unf it was also trimming off all after:
#               sequence:RefSeq:WP_014499727.1"
#                     /note="Derived by automated computational analysis using
#                     gene prediction method: Protein Homology."
# now corrected

	my $product = "NA";


	if($prefered_field_for_parsing_gene_labels == 1)
		{
		# this was the original setup. 
		# must have been a good reason for choosing these in insect mtgenomes,
		# although for plant chloroplasts, seems illogical, thus now optioned
		if($gb_feature =~ /\s+\/gene="([\w\d\s\-]+)"/){$product = $1};# 28/06/14: new regex for parsing names, previously was missing lots of 16S due to space or hyphen in each label
		if($gb_feature =~ /\s+\/note="([\w\d\s\-]+)"/){$product = $1};
		if($gb_feature =~/\s+\/product="([\w\d\s\-]+)"/){$product = $1};
		};
	if($prefered_field_for_parsing_gene_labels == 2)
		{
		if($gb_feature =~ /\s+\/note="([\w\d\s\-]+)"/){$product = $1};
		if($gb_feature =~/\s+\/product="([\w\d\s\-]+)"/){$product = $1};
		if($gb_feature =~ /\s+\/gene="([\w\d\s\-]+)"/){$product = $1};
		};


	$product = uc($product);$product =~ s/[^A-Z0-9]//g;
	my $standardized_product  = $product;
	if(exists($gene_synonyms{$product})){$standardized_product =$gene_synonyms{$product} };



# comparing complete chloroplast genomes,
# makes more sense to prefer gene field over product field 
# since, single gene label 'rpoB' has 3 product labels:
#  RNA polymerase beta subunit / RNA polymerase b-subunit / beta subunit of RNA polymerase



#############################################################

 # extract dna sequece of feature, from complete sequence

my $feature_start;my $feature_end; my $complement = 0; my $extracted_feature = "";
if($gb_feature =~ /complement\((\d+)\.\.(\d+)\)/)
	{#  CDS             complement(14968..16302)
	$feature_start = $1; $feature_end = $2;$complement = 1;
	}elsif($gb_feature =~ /(\d+)\.\.(\d+)/){
	$feature_start = $1; $feature_end = $2;
	}else
	{
	die "\nnot parsed feature posiiton:$gb_feature\n";
	};
if($feature_start =~ /\d/)
	{
#	print "feature start:$feature_start end:$feature_end\n";
	# perform some basic checks:
	if($feature_end <= $feature_start)
		{die "\nerror. why feature end ($feature_end) <= start ($feature_start) ?\n"};
	if($feature_end > length($current_dna)){print "\nerror. feature end position ($feature_end) is greater than sequence length (" , length($current_dna) , ")\n"; die ""};
	my $feature_length = $feature_end - $feature_start;
	# perl works from index 0, while ncbi reports from index 1, thus -1 here:
	$extracted_feature = substr $current_dna , ($feature_start-1) , ($feature_length+1);
	# im going to assume complement means reverse complement
	if($complement == 1){$extracted_feature =~ tr/ACGTacgt/TGCAtgca/; $extracted_feature = reverse($extracted_feature)};
	print FEATURE_DEBUG "$gi_number2\t$current_sp_id\t$proteinID\t$product\t$feature_start\t$feature_end\t$complement\t$extracted_feature\n";
	};


#############################################################




	if($gb_feature =~/\s+\/translation="(.+)"/s)# reads (AA) string over multi-lines
		{
		my $protein_sequence = $1;$protein_sequence =~ s/[\s\n\r]//g;
		$protein_sequence =~ s/^([A-Z]+)\".+/$1/;
		unless($protein_sequence =~ /^[A-Z]+$/)
			{die "\nerror 944: strange characters in protein sequence:($protein_sequence)\n"}

		if($store_proteinIDs{$proteinID} == 1)
			{
			print "warning, duplicate protein ID:$proteinID, not printing\n"
			}else{
			my $printID = $proteinID;
			if($id_format == 2){$printID = $gi_number2 . "_" . $proteinID}; # 1=protein id only, 2=gi then protein id
			if($id_format == 3){$printID = $gi_number2 . "_" . $sequence_count;$sequence_count++}; # 1=protein id only, 2=gi then protein id
			if($id_format == 4){$printID = $current_sp_id . "_" . $sequence_count;$sequence_count++};
			if($id_format == 5){$printID = $standardized_product;$sequence_count++};
			if($id_format == 6){$printID = $current_sp_id . "_" . $standardized_product;$sequence_count++};

			# print "gi_number2:$gi_number2 proteinID:$proteinID current_sp_id:$current_sp_id product:$product standardized_product:$standardized_product\n";

			if($organize_outfiles_by_genename == 1)
				{
				if($dont_print_duplicate_seq == 1)
					{
					unless(exists($store_seqs{$protein_sequence}))
						{
						open(OUT5, ">>$out_suffx.$standardized_product") || die "\nerror 506\n";
						print OUT5 ">$printID\n$protein_sequence\n";close OUT5;
						open(OUT55, ">>$out_suffx.$standardized_product.DNA") || die "\nerror 507\n";
						print OUT55 ">$printID\n$extracted_feature\n";close OUT55;

						};
					$store_seqs{$protein_sequence}=1;

					}else{
					open(OUT5, ">>$out_suffx.$standardized_product") || die "\nerror 506\n";
					print OUT5 ">$printID\n$protein_sequence\n";close OUT5;
						open(OUT55, ">>$out_suffx.$standardized_product.DNA") || die "\nerror 507\n";
						print OUT55 ">$printID\n$extracted_feature\n";close OUT55;

					};

				}else{

				if($dont_print_duplicate_seq == 1)
					{
					unless(exists($store_seqs{$protein_sequence}))
						{print USEROUT ">$printID\n$protein_sequence\n"};
					$store_seqs{$protein_sequence}=1;
					}else{
					print USEROUT ">$printID\n$protein_sequence\n";
					
					}
				};

			};
		$store_proteinIDs{$proteinID} = 1;


	#	$concatenated_protein_seqeucnes_for_this_entry .= $protein_sequence;
		$number_of_protein_features_parsed++;
		$number_of_features_found++;
	#	print "gb_feature WITH translation:$gb_feature\n";
		}else{
		$features_without_translation++;
	#	print "gb_feature WITHOUT translation:$gb_feature\n";

		};

	};

#print "number_of_features_found:$number_of_features_found
#features_without_translation:$features_without_translation\n";

if($number_of_features_found == 0){print "warning, no PROT features found in this entry\n"};


return($concatenated_protein_seqeucnes_for_this_entry);

};


#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


##########################################################################################################################################
#
#
#
#
#
##########################################################################################################################################





sub read_gene_synonyms
{

print "\nsub read_gene_synonyms\n";

# what to do with 'putative' , or 'similar to' , and misspellings: CYTOCHOMEOXIDSEII
# cad is probably short for cadherin, but what is relation to e-carherin, cadherin-1, cdh1?


%gene_synonyms = (
	'12S' 			=> '12S',
	'12SRDNA' 		=> '12S',
	'SRRNA' 		=> '12S',
'12SRIBOSMALRNA' 		=> '12S',
'SMALLSUBUNITRIBOSOMALRNA' 		=> '12S',

	'16S' 			=> '16S',
	'16SRDNA' 		=> '16S',
	'16SRIBOSOMALRNA' 	=> '16S',
	'16SRRNA' 		=> '16S',
'16SRIBOSMALRNA' 		=> '16S',

	'LRRNA' 		=> '16S',
	'LARGESUBUNITRIBOSOMALRNA'=> '16S',

	'18SRRNA' 		=> '18S',
	'18SRRNAGENE' 		=> '18S',

	'28S' 			=> '28S',
	'28SRDNA' 		=> '28S',
	'28SRIBOSOMALRNA' 	=> '28S',
	'28SRIBOSOMALRNAGENE' 	=> '28S',
	'28SRRNA' 		=> '28S',
	'28SRRNAGENE' 		=> '28S',

	'58SRRNA' 		=> '58SRRNA',
	'58SRRNAGENE' 		=> '58SRRNA',
	'5SRRNA' 		=> '58SRRNA',

	'CO1' 			=> 'CO1',
	'COI' 			=> 'CO1',
	'COX1' 			=> 'CO1',
	'COXI' 			=> 'CO1',
	'CYTOCHROMEOXIDASEI' 	=> 'CO1',

	'CO2' 			=> 'CO2',
	'COII' 			=> 'CO2',
	'COX2' 			=> 'CO2',
	'COXII' 		=> 'CO2',

	'CO3' 			=> 'CO3',
	'COIII' 		=> 'CO3',
	'COX3' 			=> 'CO3',
	'COXIII' 		=> 'CO3',

	'COB' 			=> 'CYTB',
	'CYTB' 			=> 'CYTB',
	'CYTOCHROMEB' 		=> 'CYTB',

	'CYTC' 			=> 'CYTC',
	'CYTOCHROMEC' 		=> 'CYTC',

	'EF1' 			=> 'EF1',
	'EF1A' 			=> 'EF1',
	'EF1AF1' 		=> 'EF1',
	'EF1AF2' 		=> 'EF1',
	'EF1ALPHA' 		=> 'EF1',
	'EF1ALPHA100E' 		=> 'EF1',
	'EF1ALPHAF2' 		=> 'EF1',
	'ELONGATIONFACTOR1ALPHA'=> 'EF1',
	'ELONGATIONFACTOR2LIKE' => 'EF1',

	'NAD1' 			=> 'NAD1',
	'ND1' 			=> 'NAD1',
	'NADH1' 		=> 'NAD1',
	'NDI' 			=> 'NAD1',

	'NAD2' 			=> 'NAD2',
	'ND2' 			=> 'NAD2',
	'NADH2' 		=> 'NAD2',

	'NAD3' 			=> 'NAD3',
	'ND3' 			=> 'NAD3',
	'NADH3' 		=> 'NAD3',

	'NAD4' 			=> 'NAD4',
	'NAD4L' 		=> 'NAD4',
	'ND4' 			=> 'NAD4',
	'ND4L' 			=> 'NAD4',
	'NADH4' 		=> 'NAD4',

	'NAD5' 			=> 'NAD5',
	'ND5' 			=> 'NAD5',
	'NADH5' 		=> 'NAD5',

	'NAD6' 			=> 'NAD6',
	'ND6' 			=> 'NAD6',
	'NADH6' 		=> 'NAD6',

# later additions:

'CYTOCHROMEOXIDASESUBUNIT1'  	=> 'CO1',
'CYTOCHROMEOXIDASESUBUNIT2'  	=> 'CO2',
'CYTOCHROMEOXIDASESUBUNIT3'  	=> 'CO3',

'CYTOCHROMEOXIDASESUBUNITI'  	=> 'CO1',
'CYTOCHROMEOXIDASESUBUNITII'  	=> 'CO2',
'CYTOCHROMEOXIDASESUBUNITIII'  	=> 'CO3',

'CYTOCHROMECOXIDASESUBUNITI'  	=> 'CO1',
'CYTOCHROMECOXIDASESUBUNITII'  	=> 'CO2',
'CYTOCHROMECOXIDASESUBUNITIII' 	=> 'CO3',

'CYTOCHROMECOXIDASESUBUNIT3' 	=> 'CO3',

'SIMILARTOCYTOCHROMEOXIDASESUBUNIT1'  	=> 'CO1',
'CYTOCHROMECOXIDASESUBUNIT1'  	=> 'CO1',
'SIMILARTOCYTOCHROMEOXIDASESUBUNITI'  	=> 'CO1',
'CYTOCHROMECOXIDASE1'  	=> 'CO1',
'CYTOCHROMECOXIDASEI'  	=> 'CO1',
'CYTOCHROMEOXIDASECSUBUNITI'  	=> 'CO1',
'CYTOCHROMEOXIDASE1'  	=> 'CO1',
'CYTOCHROMECOXIDASEISUBUNIT'  	=> 'CO1',

'CYTOCHROMECOXIDASESUBUNIT2'  	=> 'CO2',
'CYTOCHROMECOXIDASEII'  	=> 'CO2',

'NADHDEHYDROGENASE1'  => 'NAD1',
'NADHDEHYDROGENASE2'  => 'NAD2',
'NADHDEHYDROGENASE3'  => 'NAD3',
'NADHDEHYDROGENASE4'  => 'NAD4',
'NADHDEHYDROGENASE5'  => 'NAD5',

'NADHDEHYDROGENASESUBUNIT1'=> 'NAD1',
'NADHDEHYDROGENASESUBUNIT2'=> 'NAD2',
'NADHDEHYDROGENASESUBUNIT3'=> 'NAD3',
'NADHDEHYDROGENASESUBUNIT4'=> 'NAD4',
'NADHDEHYDROGENASESUBUNIT5'=> 'NAD5',

'NADHDEHYDROGENASESUBUNIT5ND5'=> 'NAD5',
'NADHDEHYDROGENASESUBUNITI'=> 'NAD1',

'NADHDEHYDROGENASESUBUNIT6'=> 'NAD6',

# nad4 and nad4L seem to be two seperate loci, adjacent
'NADHDEHYDROGENASESUBUNIT4L'=> 'NAD4L',

'NADHSUBUNIT1'=> 'NAD1',
'NADHSUBUNIT2'=> 'NAD2',
'NADHSUBUNIT3'=> 'NAD3',
'NADHSUBUNIT4'=> 'NAD4',
'NADHSUBUNIT4L'=> 'NAD4L',
'NADHSUBUNIT5'=> 'NAD5',
'NADHSUBUNIT6'=> 'NAD6',

'NADHDEHYDROGENSESUBUNIT1'=> 'NAD1',
'NADHDEHYDROGENSESUBUNIT2'=> 'NAD2',
'NADHDEHYDROGENSESUBUNIT3'=> 'NAD3',
'NADHDEHYDROGENSESUBUNIT4'=> 'NAD4',
'NADHDEHYDROGENSESUBUNIT4L'=> 'NAD4L',
'NADHDEHYDROGENSESUBUNIT5'=> 'NAD5',
'NADHDEHYDROGENSESUBUNIT6'=> 'NAD6',

'12SRIBOSOMALRNA' => '12S',

'CYTOCHROMEOXIDASEII'  	=> 'CO2',
'CYTOCHROMEOXIDASEIII' 	=> 'CO3',

'WINGLESS'=>'WNT',
'WINGLESSPROTEIN'=>'WNT',

'CAD'=>'CAD',
'CADHERIN'=>'CAD',

'28SLARGESUBUNITRIBOSOMALRNA' => '28S',
'12SSMALLSUBUNITRIBOSOMALRNA'=> '12S',
'16SLARGESUBUNITRIBOSOMALRNA'=> '16S',
'18SRIBOSOMALRNA' => '18S',
'18SRIBOSOMALRNAGENE'=> '18S',
'18SSMALLSUBUNITRIBOSOMALRNA'=> '18S',

'INTERNALTRANSCRIBEDSPACER1ITS1' => 'ITS1',

'INTERNALTRANSCRIBEDSPACER2ITS2' => 'ITS2',

'HISTONEH3' => 'H3',
'HISTONE3'=> 'H3',

'TRANSLATIONELONGATIONFACTOR1ALPHA'=> 'EF1',
'ELONGATIONFACTOR1A'=> 'EF1',

'ATPASE6' 	=> 'ATP6',
'ATPASESUBUNIT6' 	=> 'ATP6',
'ATPSYNTHASEFOSUBUNIT6' => 'ATP6',
'ATPSYNTHASEF0SUBUNIT6' => 'ATP6',
'ATPSYNTHASESUBUNIT6' => 'ATP6',

'ATPASESUBUNIT8' 	=> 'ATP8',
'ATPSYNTHASEF0SUBUNIT8' => 'ATP8',
'ATPSYNTHASEFOSUBUNIT8' => 'ATP8',

'CADHERINLIKEPROTEIN'=>'CAD',



# above are mostly insect mitogenome, some common nuclear genes
# below, chloroplast

# rpoB	RNA polymerase beta subunit / beta subunit of RNA polymerase / RNA polymerase b-subunit







	);


my @gene_synonym_keys = keys %gene_synonyms;@gene_synonym_keys = sort @gene_synonym_keys;

foreach my $key(@gene_synonym_keys)
{
#my $belongs_to = $gene_synonyms{$key};
#print "$key standardized to $belongs_to\n";

}



};


##############################################################################################################
#
#
#
#
#
##############################################################################################################



