#!/usr/bin/perl



###################################################################################################################################
#
#
#
# 	create_fasta_database_from_genbank_flatfiles.pl, 
#	Perl script for parsing flatfile format genbank files for certain fields, and writing in fasta format
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
# 	
# 	
#
# 	sep2009 option to ignore the list of daily releases. these are less important if you run this script infrequently
# 	oct2009 bugfix. all field variables reset after // is reached in flatfile. 
# 		bugfix: would crash if current_taxa string was found in the title but not the lineage (eg a coevol study)
# 	feb2011: simplified flatfile parsing.
# 	jul2011: take a conservative approach and do not use dna seqs with non-standard chars. print list of species which are not used.
# 	jul2012: some modification
# 	aug2013: allows multiple genbank divisions
# 	oct2013: quick hack so only specified countries are printed, see line ####unless($current_country =~ /
# 	jan2014: minor change, prints file progress
# 	9mar2014: reads the division. since you may not want est sequences, they seem to be innapropriate for species clustering, more variable than dna
# 	18jun2014: parse gi number, required by SAP species identification tool
# 	25aug2014: start using gi number instead of accession. it has more consistant format than accession 
#		(e.g. latter sometimes uses internal underscore)
# 	01oct2014: user can read in a file of accessions, only these will be output. run like this:
# 		perl create_fas... -accessions_list file.fas
#	30dec2014: parses protein sequence, from 'translation' in each feature. 
#		although hardly any translations, just a few COIs. seems not worth it.
#	30sep2015: further work on protein parsing, although might make standalone script for that
# 	13dec2015: memory efficient db reading, some of the daily releases are v big 
# 	12feb2017: switch back to accession since ncbi discontinued GI, more awkward to parse AC but whateva.
# 			improved geographic parsing, incl a number of regex's for getting co-ordinates.
# 	15may2017: added specimen voucher to metadata parsed, might need it for matching individuals
# 	01Jun2018: simple check conducted for truncated files, gives warning if found
# 	
# 	
# 	
# 	
# 	
# 	
# 	
# 	
# 	
# 	note, a gb inv file from 1 release is not the same as a gb inv from a different release. so download all for each update.
#
#
#
#
# 	wget ftp://ftp.ncbi.nih.gov/genbank/gbpln1.seq.gz --user=XXXXXX --password=XXXXXX
#
#
#
#
#
#
#



my $arg_string  = join ' ', @ARGV;


# OPTIONS:


# list of all divisions and approx number of files: 
# gbpln 1..60; gbbct 1..94; gbenv 1..57; gbinv 1..33; gbmam 1..8; gbpri 1..45; gbrod 1..30; gbvrt 1..28; gbuna 1

# put division names here, 1 or more:

my @genbank_divisions  		= 
	(

# "tempfile"
	"gbinv"
	#"gbpln","gbbct","gbenv",
#	"gbmam","gbpri","gbrod"#,"gbvrt" #,"gbuna"
#"sequence"
# "gb"#, "gbmam", #"gbpln", "gbpri", "gbrod", "gbvrt"
# "nc"#,"nc02","nc03"


	);



my $gene_specific_analysis 	= 0;	# if you have one particular gene you are interested in, put $gene_specific_analysis =1;
					# and the name of the gene in the $product_of_interest variable.
					# whereas if you are analysing many genes, put $gene_specific_analysis =0;
					# note: all entries are still printed to outfile, 
					# its just entries can be labelled either you required gene name or OtherGene

#my @product_of_interests 	= ("psbA-trnH", "trnH-psbA");
my @product_of_interests 	= ("16S");



		# 
my $print_genome_taxa 		= 0;	# 0=dont print genome taxa, 1=do print them
					# current list of these taxa, you can add more (within the function, not here)
					# Arabidopsis thaliana|Glycine max|Medicago truncatula|Oryza sativa|Populus trichocarpa|Sorghum bicolor|Vitis vinifera|Zea mays
					# Acyrthosiphon pisum|Aedes aegypti|Anopheles gambiae|Apis mellifera|Bombyx mori|Culex quinquefasciatus|Drosophila|Drosophila melanogaster|Drosophila pseudoobscura|Pediculus humanus|Tribolium castaneum
$upper_entry_length 		= 32000; # note, insect mtgenomes usually 16,000bp

$limit_taxon = 1;
$limit_taxon_name = "Hexapoda"; # Hexapoda , Insecta, Apoidea

# $parse_protein 			= 0;


my $verbose			= 0;


my $date = localtime time;$month = $date;
$month=~ s/\w+\s+(\w+)\s+\d+\s+\S+\s+(\d+)$/$1$2/;
#$month= "Apr2013";

# my $output_fasta_db_name	= "invD";


#####################################################################################################################################

$search_for_accessions = 0;
%accessions_to_look_for;

#####################################
read_command_arguments($arg_string);#
#####################################

#my $accession_keyfile_ID = $zipped_genbank_flatfiles[0];$accession_keyfile_ID =~ s/^([a-zA-Z]+).+/$1/;
my $accession_keyfile_ID = $output_fasta_db_name;



#### these old settings can be ignored:




# Hymenoptera Taxonomy ID: 7399

$download_genbank_flatfiles 		= 0; 		# 0==dont download flatfiles, 1==download. the flatfiles are updated every 3months or so
#$parse_genbank_flatfiles		= 1;
$download_daily_flatfiles		= 0;		# all taxa (not just inv). updated every day since last full release
$download_genbank_taxonomy_database 	= 0;
$ignore_daily_release			= 1;		# this differs from $download_daily_flatfiles option (which assumes these have already been downloaded, they are still parsed).

##################################################


my @zipped_genbank_flatfiles;

# scan local directory and get all file names for the divisions specified by the user:

if($user_infile =~ /[\w\d]/)
	{
	push @zipped_genbank_flatfiles , $user_infile;
	}else{
	my $file_list = "";
	foreach my $division(@genbank_divisions)
		{$file_list .= `ls $division*`}
	@zipped_genbank_flatfiles = split /\n/ , $file_list;

	#@zipped_genbank_flatfiles = ("proteobacteria_orders.gb");

	};

print "zipped_genbank_flatfiles:@zipped_genbank_flatfiles\n";


if($parse_protein 			== 1)
{
	open(PROTEIN_OUT , ">prot") || die "\n836\n";
	close PROTEIN_OUT;
};

my %sp_rm = ();
$countnumberprinted=0;
$counttotal=0;
$inv_number;
my $number_of_this_gene 	= 0;
$zipped_infileTEST ="";



#if($download_genbank_taxonomy_database==1)
#	{
#	system("rm taxdump.tar.gz");
#	$command = "wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz -O taxdump.tar.gz";print "command:$command\n";system($command);
#	system("rm *.dmp");
#	system("tar xvzf taxdump.tar.gz");
#	}else{print "NOT downloading taxonomy database\n"};


# tar xvzf taxdump into ~/hip-db/hip-db_scripts/ folder, delete all but names.dmp and nodes.dmp.
# copy coleoptera_fasta to ~/hip-db/ folder
# ftp://ftp.ncbi.nih.gov/blast/executables/LATEST/ for ncbi tools

# improvements:should look for holmetabola as well as endop, just in case




$outfile1_name = "database_fasta_" . $month;





############################
parse_genbank_flatfiles2();#
############################




exit;





##################################################################################

print "downloading list of daily update flatfiles\n" , $month , "\n";
$command = "rm flatfilelist_" . $month;print "command:$command\n";system($command);
$command = "wget ftp://ftp.ncbi.nih.gov/genbank/daily-nc/ -O " . "flatfilelist_" . $month;print "command:$command\n";system($command);



unless($ignore_daily_release==1)
	{


	open(IN, "flatfilelist_$month") || die "cant open flatfilelist_$month\n";
	while ($line= <IN>)
		{#/genbank/daily-nc/nc0219.flat.gz"
		if($line=~ /\/genbank\/daily\-nc\/(nc\d+\.)flat\.gz\"/)
			{

			if($download_daily_flatfiles == 1)
				{
				$command = "rm $1$month.gz";print "command:$command\n";system($command);
				$command = "rm $1$month";print "command:$command\n";system($command);
				$command = "wget ftp://ftp.ncbi.nih.gov/genbank/daily-nc/" . $1 . "flat.gz -O " . $1 . $month . ".gz";	print "command:$command\n";	system($command);
				$command = "gunzip " . $1 . $month . ".gz";print "command:$command\n";system($command);
				};



			$infile = $1 . $month;
	#		parse_genbank_flatfile();
			system("rm $infile");
			
			# $infile_names[$count_inv_files] = $1 . $month;$count_inv_files++;


			}
		}
	close IN;

	}

##################################################################################


print "downloading list of release flatfiles\n" , $month , "\n";
$command = "rm releaseflatfilelist_" . $month;print "command:$command\n";system($command);
#$command = "wget ftp://ftp.ncbi.nih.gov/genbank/ -O " . "flatfilelist_" . $month;print "command:$command\n";system($command);
$command = "wget ftp://ftp.ncbi.nlm.nih.gov/genbank/ -O " . "releaseflatfilelist_" . $month;print "command:$command\n";system($command);

$debugging_counter=0;
$printed_counter=0;
$entry_counter=0;
$title_counter = 1;



$count_inv_files=0;
open(IN, "releaseflatfilelist_$month") || die "cant open releaseflatfilelist_$month\n";
while ($line= <IN>)
	{
	if($line=~ /(gbinv\d+\.)/)
		{
		print $line;

		if($download_genbank_flatfiles == 1)
			{
			$command = "rm $1$month.gz";print "command:$command\n";system($command);
			$command = "rm $1$month";print "command:$command\n";system($command);
			$command = "wget ftp://ftp.ncbi.nih.gov/genbank/" . $1 . "seq.gz -O " . $1 . $month . ".gz";	print "command:$command\n";	system($command);
	#		$command = "gunzip " . $1 . $month . ".gz";print "command:$command\n";system($command);
			};

		$infile = $1 . $month;
#		parse_genbank_flatfile();
#		system("rm $infile");


		$infile_names[$count_inv_files] = $1 . $month;$count_inv_files++;

		};
	};
close(IN);
$command = "rm flatfilelist_" . $month;print "command:$command\n";system($command);

#close(OUT);
#close(OUT2);
#close(OUT3);

print "\nEND OF SCRIPT\n";
exit;


















##############################################################################################




sub parse_genbank_flatfiles2
{
open(LOG2, ">>DodgyDnaSeqsFound") or die "cantopen\n";
open(MAIN_OUT, ">$output_fasta_db_name") or die "cantopen\n";
open(OUTKEY, ">metadata.$month.$accession_keyfile_ID") || die "\nerror 404\n";
print OUTKEY "accession\ttaxid\torder\ttaxstring\tspecies_id\tproduct\tcountry\tcoordinates\n";

close(OUTKEY);

print LOG2 "\n\n\nRUNNING:$date\n";


for my $zipped_infile(@zipped_genbank_flatfiles)
	{
	$countfiles++;
	#$zipped_infile = "gbpln51.seq.gz";
	#my $command = "tar -xvzf $zipped_infile";

	my $command = "gunzip $zipped_infile";
	print "$countfiles of $#zipped_genbank_flatfiles command:$command\n";system($command);

	$zipped_infile =~ s/\.gz//;
	$zipped_infileTEST = $zipped_infile;

	####################################
	parse_genbank_flatfileNEW($zipped_infile);#
	####################################

	my $command = "gzip $zipped_infile";print "command:$command\n";system($command);

	print "\nfrom all files, $countnumberprinted were printed out of $counttotal ..... ($num_printed2)\n";
	print "n_protein_features_parsed:$number_of_protein_features_parsed n_entries_with_protein:$number_of_entries_with_protein\n";
	if($gene_specific_analysis ==1){print "\n\nnumber of @product_of_interests is $number_of_this_gene\n"}

	}



close(LOG2);
close(MAIN_OUT);

my @speciesremoved = keys %sp_rm;
print scalar @speciesremoved , " removed\n@speciesremoved\n";
print "\n$countnumberprinted were printed out of $counttotal\n";

if($truncated_files =~ /[\w\d]/)
	{
	print "warning, following files are not terminated as expected, check they are downloaded completly:\n";
	print "$truncated_files\n";
	};


}




###########################################################################

sub parse_genbank_flatfileNEW
{
my $inputfile = shift;
print "\n\nsub parse_genbank_flatfileNEW\n\n";




open ( IN2 , $inputfile ) || die "cant open infile:$inputfile\n";
print "opened $inputfile\n";

#my $fileasstring = "";
#while(my $line = <IN2>)
#	{$fileasstring .= $line}
#my @file_as_array = split( /\n\/\/\n/ , $fileasstring );	# problem here, crashed with a big file, 
							# not sure if its just the size of the file or some feature of.


$current_file_entry_index = 0;							# perl error message was 'split loop'

my $fasta_entry = "";my $last_line;
while(my $fasta_line= <IN2>)
	{
	if($fasta_line =~ /^\/\/\n*$/)
		{
		unless(length($fasta_entry)<=2)
			{
			$current_file_entry_index++;

			unless($print_genome_taxa == 0 && $current_length >= $upper_entry_length)
				{
				########################################
				parse_this_entry($fasta_entry, $inputfile);#
				########################################
				}
			}
		$fasta_entry = "";
		}else{
		my $current_length = length( $fasta_entry );
		unless($current_length >= 1000000000)
			{
			$fasta_entry .= $fasta_line;
			};
		}
	$last_line = $fasta_line;
	};

close(IN2);


# print "last_line:($last_line)\n"; 

if($last_line =~ /^\/\//)# file terminated with expected delimiter
	{
#	print "file looks ok\n";	
	}else{
	print "warning, file not terminated as expected, possibly not completly dowloaded, check file:$inputfile\n";
	$truncated_files .= "$inputfile ";
	};





foreach my $entry(@file_as_array)
	{
	#print "\n$inputfile. from all files, $countnumberprinted were printed out of $counttotal\n";
	#if ($entry =~ /Z93710/i){print $entry}
#if($entry =~ /KC880663|JQ024878|HQ266313/){
	#########################
#########	parse_this_entry($entry);#	
	#########################
#die};
	}



}





###########################################################################



sub parse_this_entry
{
my $current_entry = $_[0]; my $current_filename = $_[1]; #print "current_entry:$current_entry\n";die;
my $current_entry_copy = $current_entry;
$counttotal++;

#if($zipped_infileTEST =~ /inv4/){print "\n\ncurrent_entry:$current_entry"}

my $year = "unknown_year";
if($current_entry =~ /JOURNAL\s.+\D([12][78901][0-9][0-9])\D.*\n/)
	{
	$year = $1;$year = "P$year"; # publication year
	}elsif($current_entry =~ /LOCUS\s+.+\s\d+\-[A-Z]+\-([12][78901][0-9][0-9])\s*\n/)
	{
	#print "WARNING:cannot find year in JOURNAL line of this entry: $current_entry";
	$year = $1;$year = "M$year"; # modification year
	}else{print  "WARNING:cannot find year in LOCUS line of this entry: $current_entry"}
# LOCUS       SCU49845     5028 bp    DNA             PLN       21-JUN-1999


#JOURNAL   Genes Dev. 10 (7), 777-793 (1996)
#JOURNAL   Submitted (22-FEB-1996)



# division (since daily releases can contain est/rna etc)

my $inv_division = 0;

if($current_entry =~ /(LOCUS\s+.+INV\s\S+)\n/)
	{
	my $locus_line = $1;
	#print "locus_line:$locus_line\n";
	$inv_division=1;
#LOCUS       HX672325                 687 bp    mRNA    linear   EST 09-FEB-2014
	}




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
		print "\n\nERROR, no accession\n\n$current_entry\n"
		}
	}

#VERSION     U49845.1  GI:1293613
my $gi_number = "UnknownGI";
if($current_entry =~ /VERSION\s+\S+\s+GI\:(\d+)\n/)
	{
	# single values for accession found
	$gi_number = $1
	};

#print "\n\ngi_number:$gi_number\ncurrent_entry:$current_entry\n";

#######		HEXAPODS

my $in_hex = 0;# && $current_entry =~ /18S|small\s*subunit/i
#if($current_entry =~ /hexapoda/i){$in_hex=1}
if($current_entry =~ / $limit_taxon_name\;/i){$in_hex=1};


#            Eukaryota; Metazoa; Arthropoda; Hexapoda; Insecta; Pterygota;
#            Neoptera; Paraneoptera; Hemiptera; Sternorrhyncha; Aphidiformes;
#            Aphidoidea; Aphididae; Aphidini; Aphis; Aphis.

if($current_entry =~ /Coleoptera\;/){$number_of_coleoptera++}
if($current_entry =~ /Tribolium castaneum/){$number_of_Tribolium++}

#######		ORDER

my $order = "NA";
while($current_entry =~ /([A-Z][a-z]+tera)\;/){$order = $1;$current_entry =~ s/[A-Z][a-z]+tera\;//}
#print "$order\n";




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
		print "\n\nERROR, no tax id\n\n$current_entry\n\n";
		
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
my $discard_entry = 0;
my $bold = 0;




#######		ORGANISM

	if($current_entry =~/ORGANISM\s*(.+)\n/)
		{
		$current_organism = $1;$current_taxstring = $current_organism;$current_taxstring =~ s/^\s//g;$current_taxstring =~ s/\s$//g;$current_taxstring =~ s/\s/_/g;
#print "current_taxstring:$current_taxstring\n";
	#	if($current_organism=~/Drosophila melanogaster|Bombyx mori|Tribolium castaneum|Apis mellifera|Anopheles gambiae/i)
	#	if($current_organism =~ /Acyrthosiphon pisum|Aedes aegypti|Anopheles gambiae|Apis mellifera|Bombyx mori|Culex quinquefasciatus|Drosophila|Drosophila melanogaster|Drosophila pseudoobscura|Pediculus humanus|Tribolium castaneum/i)

		# NOTE, these are ignored later if you have specfied to print genome taxa
		# so no need to comment out
		if($current_organism =~ /Arabidopsis thaliana|Glycine max|Medicago truncatula|Oryza sativa|Populus trichocarpa|Sorghum bicolor|Vitis vinifera|Zea mays/i || 
			$current_organism =~ /Acyrthosiphon pisum|Aedes aegypti|Anopheles gambiae|Apis mellifera|Bombyx mori|Culex quinquefasciatus|Drosophila|Drosophila melanogaster|Drosophila pseudoobscura|Pediculus humanus|Tribolium castaneum/i
			){$discard_entry=1};
		if($current_organism=~/Homo sapiens/i)
			{$discard_entry=1};

		if($current_organism =~ /BOLD\:/){$bold = 1
			#;$current_organism =~ s/BOLD\:/BOLD/;
			}

		}

	if($current_organism =~/BOLD/ && $bold == 0)
		{
	#	print "\nerror 602 BOLD:$current_entry\nerror 602...\n"
	# not currently used
		}


#######		PRODUCT

my $current_product = "NA";
	if($current_entry =~/\/product..(.+)\"/)
		{
		$current_product = $1;$current_product =~ s/[^a-zA-Z0-9]//g;
		};




#######		SPECIES ID

my $current_sp_id = "NA";

	if($current_entry =~/\/organism\=\"([^\"]+)\"/)
		{#/organism="Catenula sp."
		$current_sp_id = $1;$current_sp_id =~ s/[\n\r\.]//g;#$current_sp_id =~ s/^([a-zA-Z\s]+)[\:\,].+/$1/g;
		};



my $specimen_voucher = "NA";
if($current_entry =~/\/specimen_voucher\=\"([^\"]+)\"/)
	{
	$specimen_voucher = $1;
	};




#######		COUNTRY

my $current_country = "NA";
my $current_region = "NA";


	if($current_entry =~/\/country..([^\"]+)\"/)
		{
		$current_country = $1;$current_country =~ s/[\n\r]//g;
		if($current_country =~ s/^([a-zA-Z\s]+)[\:\,](.+)/$1/)
			{$current_region = $2; $current_region =~ s/  +/ /g};
#	                     /country="Brazil: Sao Paulo state, Pariquera Acu
#	                     municipality, Pariquera Mirim district"
#	                     /lat_lon="24.71667 S 47.86667 W"
		};


#######		COORDINATES

# unfortunately, coordinates are really inconsistently written by uploaders. 
# i try to have specific regexs to avoid picking up non-coordinates as much as possible, 
# hence many instead of small number of generalist ones.

my $current_coordinates = "NA";

# /note="Geographic coordinates: 28:20N/95:32W;
# /lat_lon="09deg 29'N, 99deg 57'E"



	if($current_entry =~ /lat_lon=\"([\d\.]+ [nsNS] [\d\.]+ [ewEW])\"/)#/lat_lon="35 N 83 W"
			{
			$current_coordinates = $1;
			}elsif($current_entry =~ /lat_lon=\"([^\"]+[nsNS][^\"]+[ewEW])\"/)
			{
			$current_coordinates = $1;
			}elsif($current_entry =~/(\d{1,2}[\.\:]\d+\s*[nsNS][\,\/]*\s*\d{1,3}[\.\:]\d+\s*[ewEW]+)/)
			{
			$current_coordinates = $1;
			}elsif($current_entry =~/lat_lon=\"(\d+ degrees [^\"]+)\"/)
			{			# lat_lon="20 degrees 06' 00', 100 degrees 36' 00'"
			$current_coordinates = $1;
			}elsif($current_entry =~/lat_lon=\"[\d\.]+[nsNSewEW] [\d\.]+[nsNSewEW]\"/)
			{			# /lat_lon="14.63888E 77.68054N"
			$current_coordinates = $1;
			}elsif($current_entry =~/country=\"[^\"]+ ([\d\.\-]+ [nsNSewEW]\, [\d\.\-]+ [nsNSewEW])\"/) 
			{			# /country="Peru:Cuzco, Paucartambo, Consuelo, 15.9 km SW Pilcopata, 1000 m, -13.01417 S, -71.29511 W"
			$current_coordinates = $1;
			}elsif($current_entry =~/lat_lon=\"([nsNSewEW][^\"]+\d+)\"/) 
			{			# /lat_lon="N27.151667, E120.47167"
			$current_coordinates = $1;
			}elsif($current_entry =~ /lat_lon=\"(\d+[0-9\.\:\-]+\d+)\"/) 
			{			# /lat_lon="23:11.068-51:10.206"
			$current_coordinates = $1;
			}elsif($current_entry =~ /\/note=\"[\w\s]+co\-*ordinates(.[^\"]+[nsNSewEW][\;\,][^\"]+[nsNSewEW])\"/) 
			{			# /note="approximate coordinates: 22 46'12 S; 43 18'36 W"
			$current_coordinates = $1;
			}elsif($current_entry =~ /lat_lon\s+[\:]+\s*([nsNSewEW\s\.0-9]{13,19}\S)/) 
			{			# lat_lon             :: S 20.917 E167.08 3
			$current_coordinates = $1;
			}elsif($current_entry =~ /isolation_source\=\"[a-zA-Z\,\s]+[\;]([\s\d\,\.sSnN]{10,21}[eEWw])\"/) # /isolation_source="Ancud, Chiloe, Chile; 41 53 S, 73.50W"
			{
			$current_coordinates = $1;
			}elsif($current_entry =~ /\s(\d{1,3}\.\d{1,3}\s[sSNn]\s\/\s\d{1,3}\.\d{1,3}\s[WwEe])\.\"/) # # /isolation_source="Southwest Pacific, station COOK14MV-38. 15.25 S / 172.14 W."
			{ 
			$current_coordinates = $1;
			}elsif($current_entry =~ /lat_lon\=\"([^\"]{24,44})\"/) # # /lat_lon="Lat.: 48o 33.282'N; Long.: 122o 56.525'W'"
			{ 
			$current_coordinates = $1;
		




      				#   lat_lon field, or sometimes has it in the country field a la:   /country="Antarctica: 66.01S 146.24E, Mertz Glacier
#		unless($current_entry =~ /\/isolation_source\=\"|\/lat_lon\=|\/country..[a-zA-Z\:\, ]+\d+\.\d+\s*[nsNS]\,*\s*\d+\.\d+\s*[ewEW]+/)
#			{
#		#	print LOG2 "$current_entry \nabove appears to have coordinates, but not in the standard field\n\n";
#		#	die "$current_entry \nabove appears to have coordinates, but not in the standard field\n\n" 
#			}

			}elsif($current_entry =~ /(\d{1,2}\s*[nsNS]\W+\d{1,3}\.\d+\s*[ewEW])/)
			{
			my $test = $1;
			die "entry:$current_entry\ntest:$test\ndifculy psring coordinatres from above...\n";
			}elsif($current_entry =~ /lat_lon|co\-*ordinate.+\d\s*[nsNS].+\d+\s*[ewEW]/)
				{die "\n$current_entry\npretty dure entry has coordinates not paersed\n"};


# check specific fields for gene names



if ($gene_specific_analysis ==1 && length ($current_entry) <= 100000)	# for species dense / gene specific analyses, genome entries are not of interest
	{
	#print "\n\n\nentry:$current_entry\n";
	my $string_to_find_genes = "";
	#if($string =~ /blah1.+blah2/s)
	if(   $current_entry =~ /DEFINITION(.{5,800})ACCESSION/s)
		{$string_to_find_genes .= $1
		}else{
		print "\n$current_entry\n error 550. cant read definition line in entry above. excessive lines (>800 chars) are not read \n"
		}
	while($current_entry =~ /gene=\"([^\"]{1,50})\"/)
		{$string_to_find_genes .= $1;$current_entry =~ s/gene=\"([^\"]{1,50})\"//;
		if(length($string_to_find_genes) >= 100){break }
		}
	while($current_entry =~/\/product..(.{1,50})\"/){$string_to_find_genes .= $1;$current_entry =~ s/\/product..(.{1,50})\"//;
		if(length($string_to_find_genes) >= 100){break }
		}
	while($current_entry =~/\/note\"(.{1,50})\"/){$string_to_find_genes .= $1;$current_entry =~ s/\/note\"(.{1,50})\"//;
		if(length($string_to_find_genes) >= 100){break }
		}

	#print "\nstring_to_find_genes:$string_to_find_genes\n\n";

	my $found_product = 0;
	foreach my $product_of_interest(@product_of_interests)
		{
		if($string_to_find_genes =~ /$product_of_interest/i)
			{$current_product = $product_of_interest;$found_product=1};
		}
	if($found_product ==1)
		{$number_of_this_gene++}else{$current_product = "OtherGene"}

	#print "current_product:$current_product\n";
	}

my %current_gene_hash = ();

if ($gene_specific_analysis ==0 && length ($current_entry) <= 100000)
	{
	my $numberofgenes = 0;
	while($current_entry =~ /gene=\"([^\"]{1,50})\"/)
		{
		my $currentgene1 = $1;$numberofgenes++;
		$currentgene1 =~ s/[^a-zA-Z0-9]//g;
		$current_entry =~ s/gene=\"([^\"]{1,50})\"//;
		if($numberofgenes >= 20){break }else{$current_gene_hash{$currentgene1} = 1}
		}

	}
 
# gene feild should be more standardized.
my @currentgenelist = keys %current_gene_hash;
my $currtgenestring = "";if($currentgenelist[0] =~ /[\w\d]/){$currtgenestring = join(' ', @currentgenelist)}else{$currtgenestring = "NA"}



### PROTEIN

my $protein_seq = "";
if($parse_protein == 1)
	{		#####################################################
	$protein_seq = parse_protein_seqs_in_this_entry($current_entry);#
	};		#####################################################




#######		DNA

my $current_dna = "";

#$current_entry =~ s/^[\n\r\s.]+ORIGIN\s*[\n\r]//;
	$current_entry =~ s/\012\015?|\015\012?//g;

#                  /number=4
#ORIGIN      1 bp upstream of SalI site.
#        1 gtcgactgca ctcgccccca cgagagaaca gtatttaagg agctgcgaag gtccaagtca
   
$current_entry =~ s/(ORIGIN\s+)\d bp upstream of [a-zA-Z]{1,10} site\./$1/;


#print $current_entry;
#$line =~ s/^.+\n//;
	if($current_entry =~/^.+ORIGIN(.+)$/)
		{
		$current_dna = $1;$current_entryCOPY2 = $current_entry;
		$current_entry =~ s/^.+ORIGIN//;
		if($current_entry =~/^.+ORIGIN(.+)$/){		print"\n\nERROR, 2 * ORIGIN?\n\n$accession\n$current_entry_copy"; 	die "\n\nerror 423\n"}

		# bizarely slow:
		#while($current_dna =~ /\d+1\s([actgn]{10})/){$current_dna =~ s/\d+1\s[actg]{10}//}

		$current_dna =~ s/\d//g;$current_dna =~ s/[\s ]//g;
		if($current_dna =~ /^[atcgurykmswbdhvn]+$/)
			{
		#	$accession
		#	$current_taxid
		#	$current_dna

		#	unless($discard_entry == 1 || $bold == 1 || $in_hex == 0)
	#	unless($discard_entry == 1  || $in_hex == 0)

	if($print_genome_taxa == 1){$discard_entry = 0}





############################################################################################################


# limit  by country or taxon:

####unless($current_country =~ /china/i){$discard_entry =1}#####################################

if($limit_taxon == 1)
	{
	unless($in_hex == 1){$discard_entry = 1};#####################################
	};

if(length($current_dna)<=150){$discard_entry =1;$tooshort++;#print "too short ($tooshort)\t"
				};
if($accession =~ /UnknownAccesion/){$discard_entry =1;print "no accession\t"}
#unless($inv_division== 1){$discard_entry =1}#####################################

if($print_genome_taxa 		== 0)
	{if(length($current_dna)>=$upper_entry_length){$discard_entry =1}};


############################################################################################################



	if($discard_entry == 1 )
		{
		#print "not printing genomes\n";
		}else
				{
				$countnumberprinted++;
				if($taxid_present == 1)
					{
					# $current_taxid is the ncbi taxon number
					# parse_taxon_from_fastafile expects to see the accession followed by the ncbi tax number.


					if($search_for_accessions == 1 )
						{
						if($accessions_to_look_for{$accession} == 1)
							{
							# print fasta output:
							print MAIN_OUT ">$accession" , "_$current_taxid" , "\n$current_dna\n";
							# or .gb output:
							#print MAIN_OUT "$current_entry_copy\n" ."//" . "\n";$num_printed2++;
							}else{
							}
						}else{

						if($output_ID_format == 1)
							{print MAIN_OUT ">$accession" , "_$current_taxid" , "\n$current_dna\n"};

						if($protein_seq =~ /\w/ && $parse_protein == 1)
							{
							open(PROTEIN_OUT , ">>prot") || die "\n836\n";
							print PROTEIN_OUT ">$accession" , "_$current_taxid" , "\n$protein_seq\n";
							close PROTEIN_OUT;$number_of_entries_with_protein++;
							}else{};

						# $current_taxstring is the organism name:
						# gradualy switching to GI, easier to parse than accession
						if($output_ID_format == 2)
							{
						# you prob want this one for invert mtgenomes:
						 print MAIN_OUT ">$current_taxstring" , "_$accession" ,"\n$current_dna\n";	
							};

						
						}

					



					}else{
					print MAIN_OUT ">$accession" , "_$current_taxstring" ,   "__$current_product\n$current_dna\n";				
				#	print "\nprinting >$accession\n";
					print "\ntax id missing\n";
					}

				open(OUTKEY, ">>metadata.$month.$accession_keyfile_ID") || die "\nerror 586\n";
				print OUTKEY "$current_filename" , "_$current_file_entry_index\t$accession\t$gi_number\t$current_taxid\t$order\t$current_taxstring\t$current_sp_id\t$current_product\t$current_country\t$current_region\t$current_coordinates\t$year\t$currtgenestring\t$specimen_voucher\n";
				close(OUTKEY);


				}
			}else{
			my $current_dnaCOPY = $current_dna;
			$current_dnaCOPY =~ s/[atcgurykmswbdhvn]//g;
			print LOG2 "non standard chars in dna seq $accession file $inv_number chars:$current_dnaCOPY\n"; 
			$current_dna =~ s/[^atcgurykmswbdhvn]/n/g;
			$sp_rm{$current_organism}++;
			};

		}else{
		print"\n\nERROR, no DNA\n\n$current_entry\n"; 
		#die 
		
		}

#if($accession =~ /AB181702/){print $current_entry_copy}

#if($entry =~ /KC880663|JQ024878|HQ266313/){
#print  "$accession\t$current_taxid\t$order\t$current_taxstring\t$current_sp_id\t$current_product\t$current_country\t$current_coordinates\t$year\t$currtgenestring\n";
#die;
#};

  
if ($verbose == 1){print "$accession" , "_$current_taxid \n" };



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
     ***** create fasta database from genbank flatfiles *****
		  
\n";


# perl parse_clades.pl -treefile RAxML_bipartitions.diet_group.0.reformatted -seqfile diet_group.0.fas -keyfile key_Mar2013_Magnoliophyta
#$treefile;
#$fas_file;
#$output_filename 	= "$treefile.diet_clades";



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

	}

if($arguments =~ /-in\s+(\S+)/)
	{
	$user_infile = $1;
	}

if($arguments =~ /-out\s+(\S+)/)
	{
	$output_fasta_db_name = $1;
	}else{
 	$output_fasta_db_name	= "default_outfile_name";
	};



if($arguments =~ /-outformat\s+(\S+)/)
	{
	$output_ID_format = $1;
	}else{
	die "\nerror. please give output format (fasta IDs) in command, thusly:
-outformat 1 if you want accession_ncbiTaxnumber
-outformat 2 if you want Taxstring_GInumber
\n";
	}
#$output_ID_format = 1;# 1=$accession_$ncbiTaxnumber
#			# 2= $currentTaxstring_$gi_number for invert mt genomes

if($arguments =~ /-parse_protein/)
	{
	$parse_protein = 1;
	}else{
	$parse_protein = 0;
	};




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
my $gb_entry = shift;
#print "\nsub parse_protein_seqs_in_this_entry:$gb_entry\n";
my @split_gb_entry = split /\s+CDS\s+|\s+tRNA\s+|\s+rRNA\s+|\s+misc_feature\s+/ , $gb_entry;# 
print "split_gb_entry:$#split_gb_entry\n";

my $concatenated_protein_seqeucnes_for_this_entry = "";
my $number_of_features_found=0;
my $features_without_translation=0;
for  my $k(1 .. $#split_gb_entry)
	{
	my $gb_feature = $split_gb_entry[$k];#	print "\n\nFEATURE:($gb_feature)\n";
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


	if($gb_feature =~/\s+\/translation="(.+)"/s)# reads (AA) string over multi-lines
		{
		my $protein_sequence = $1;$protein_sequence =~ s/[\s\n\r]//g;
		$protein_sequence =~ s/^([A-Z]+)\".+/$1/;
		unless($protein_sequence =~ /^[A-Z]+$/)
			{die "\nerror 944: strange characters in protein sequence:($protein_sequence)\n"}
			#print "\tprotein_sequence:($protein_sequence)\n";
			#$store_current_proteins{$protein_sequence} = $product;$products_found++;
		$concatenated_protein_seqeucnes_for_this_entry .= $protein_sequence;
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



