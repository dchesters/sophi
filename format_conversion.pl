#!/usr/bin/perl

###################################################################################################################################
#
#
#
# 	format_conversion.pl, 
#	Perl script for converting between dna sequence formats
#
#    	Copyright (C) 2013-2017 Douglas Chesters
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
#	change log:
# 	sep2013: reads phylip in order, thus prints out in same order, better for specifying taxa in e.g. tnt
#	oct2013: added gpl. no longer called dougseq.pl .. 
#	jul2014: trims off numbers that pynast appends to fasta ID
#	02oct2014: error messg when you tell it to read phylip, but dont give it a phylip format file
#	06mar2015: last entry in file was not being cleaned up 
#	04DEC2015: bugfix write_fasta, last entry was previously ommited.
#	09SEP2016: checks input fasta file is actually an alignment (all entries have seq of same length),
#			error message if not.
#	02OCT2016: error message if duplicate IDs in phylip outfile
#	27APR2017: error message if read fasta invoked without fasta file
#	19MAY2017: option to remove bad bases, 
#		(i found a few J's in the protein seqs of hexapod mitigenomes, which breaks raxml)
#
#
#
#
#	notes:
#
#
#
#
#
#
###################################################################################################################################



$input 		= $ARGV[0];
$output 	= $ARGV[1];
$input_format 	= $ARGV[2];
$output_format 	= $ARGV[3];
#$output_r 	= $output . ".randomised.fas";
#$output_a 	= $output . ".alphabetized.fas";

$missing_data_character = "N";
$parse_taxa 	= 0;	# this should be switched off by default
$replace_id_spaces = 0;	# only works when reading fasta

$remove_bad_bases = 1;

# recode gaps. calculate lengths. remove spaces when reading other formats. lowercase / uppercase option

#@taxa_parse_list = (
#"CA1C3C4C7","CA1C3Ca4CHa5Pt6Pt7","CA1C3Ca4CTrBe6Be7","CA1D3D4A5Ag6Ag7","CA1D3D4A5Ag6Ilb7","CA1D3D4Co5Co6Co7","CA1D3D4Co5Co6Pa7","CA1D3D4D5Cy6Cy7",
#"CA1D3D4H5Bi6Lib7","CA1D3D4H5Hy6De7","CA1D3D4H5Hy6Hy7","CA1D3D4H5Hy6Ni7","CP1C2Ch3C4Ch5Do6Ca7","CP1C2Ch3C4DoDo7","CP1C2Ch3C4DoPl7","CP1C2Ch3C4Ga5Al6BlBl7",
#"CP1C2Ch3C4Ga5Al6unPs7","CP1C2Cuu3C4Ce5Ce7","CP1C2Cuu3C4Cu5Cu6Cu7","CP1C2Cuu3C4Sc5Ip7","CP1S2S3S4StSt5St7","CP1S2S3S4TaAl5At7","CP1Sc2S3S4Ap5Ap7","CP1Sc2S3S4Ru5An7",
#"CP1Sc2S3S4Sc5Hei7","CP1Sc2S3S4Sc5Onh7","CP1Sc2S3S4Se5Ma7","CA1C3Ca4CHa5Ha6","CA1C3Ca4CHa5Pl6","CA1D3D4D5Ac6","CA1D3D4H5Hyh6","CA1D3D4L5La6","CP1C2Ch3C4Ga5Al6"
#);



$screen_print = 1;

if($screen_print == 1)
	{
print "\n\n\n************************************************\n\n";

if($parse_taxa == 1)
	{
	print "\nyou have chosen to parse the taxa:\n@taxa_parse_list\n";
	};

print "\nscript\n";
print "reads:\tfasta nexus clustal phylip\n";
print "writes:\tfasta nexus phylip fasta_degapped fasta_missing_all fasta_missing_terminal nexus_concise(no name block) \n";
print "\tpaup_bootstrap paup_distance paup_parsimony\n";

 print "\nrequires no external perl modules.\n";
 print "assumes no redundent id's\n";
 print "also assumes no spaces in sequence id's (when reading nexus, fasta should not be affected)\n";
 print "for tnt you will probaly need nexus_concise format\n";
 print "removes length info from end of fasta id's (if input and output format both fasta), eg \">tax1 680 bp\" or \">tax1\/1-762\"";
 print "\nwill work across platforms (although this is generally bad practice). ";	
print "\nscript can remove gaps (fasta_degapped), change terminal gaps to missing characters (fasta_missing_terminal) ";
print "\n\tor change all gaps to missing (fasta_missing_all format)";
print "line 10 sets the missing data character\n\t(put N here if necessary). default is ?\n";
print "also produces output with randomised id's, only works for fasta output\n";
print "randomised output for control analyses only!\n";
print "also makes alphabetized fasta file although not quite perfect yet.\n";
print "phylip format has been tested for paml and raxml\n";
print "if converting to paup format, appropriate command block will also be printed to outfile ";
print ".\nalso if one of the taxa has the string: outgroup in the id, paup will be informed that this taxon is the outgroup\n";
print "when reading nexus format, ALL COMMENTS ARE REMOVED. this means everything between (and incl) square brackets.\n";
print "it seemed to have an issue with clustal conversion at one point but i havnt been able to replicate the problem\n";
print "paup does not like dashes, spaces or brackets in taxon names, unless quoted. (done automatically)\n";
	};

if(length($input)<=1 || length($output)<=1 || length($input_format)<=1 || length($output_format)<=1)
	{
print "\n\n\nYOU HAVE NOT TYPED COMMAND CORRECTLY. TRY AGAIN\n\n";
print "to run:\nperl format_conversion.pl input_file output_file input_format output_format\n";
die;
	};

if($screen_print == 1)
	{
print "\n\nyou have chosen:\ninput file:$input ($input_format format) output file:$output ($output_format format)\n\n";
	};

# globals
my $entry_counter = 0;
my $current_seq_length;
my @seqs;
my @ids;
my %seq_hash;
my %ids_r;
my $print_warning_message;

	# INPUT 

if ($input_format =~ /fasta/i)
	{
	read_fasta();
	};

if ($input_format =~ /nexus/i)
	{
	read_macclade();
	};

if ($input_format =~ /clustal/i)
	{
	read_clustal();
	};

if ($input_format =~ /phylip/i)
	{
	read_phylip();
	};

if ($input_format =~ /macclade/i)
	{
	read_macclade();
	};


	# OUTPUT

if ($output_format =~ /nexus/i)
	{
	write_nexus();
	};

if ($output_format =~ /fasta/i)
	{
	write_fasta();
	};

if ($output_format =~ /phylip/i)
	{
	write_phylip();
	};

if ($output_format =~ /paup/i)
	{
	write_paup();
	};





if($screen_print == 1)
	{
print "\nend of script. thankyou for using ... \n";
print "\n\n*******************************************\n\n";
	};

if($print_warning_message == 1)
	{print "\n\nWARNING: nexus file appears to be sequential yet ";
	print " some sequence id's are duplicate.\nmay be due to non-unique id's. this will cause probs with output file.\n"};


if($remove_bad_bases == 1 && $bad_bases_removed >= 1)
	{
	print "
warning, you have set script to remove bad bases, and $bad_bases_removed have indeed been found.
";

	};	




print "\n\nformat conversion script finished.\n";
exit;




####################################################################################
#
#
#
#
#
####################################################################################



sub read_fasta
{

if($screen_print == 1)
	{
	print "subroutine to read fasta format.\nreading input file ($input).\n";
	};

open(FASTA_IN, $input) || die "Cant open input:$input.\n";
print "\t.... opened ....\n";



	# local variables:

	my $line;
	my $current_id;
	my $current_sequence;
	my $current_name;
	my $current_seq;


	my $file_as_string = ""; my @all_lines = ();
	while($line= <FASTA_IN>){$file_as_string .= $line};
	close(FASTA_IN);

	unless($file_as_string =~ />/){die "\nyou specified read fasta, but input file is not fasta format. quitting.\n\n"};

	@all_lines = split /\012\015?|\015\012?/, $file_as_string;
	foreach $line(@all_lines)
		{
#		print "line:($line)\n";
#	while($line = <FASTA_IN>)
#		{ #cb1

		$line =~ s/\n//;
		$line =~ s/\r//;
		if ($line=~/^\s{0,2}>\s{0,2}(\S.*\S)\s*$/ )
			{
			print ""; #> 10CDELA4arenicola
			$line =~ s/(\S+)\s\d{1,4}\s+bp\s*$/$1/; # remove 688 bp from end of id if present
			$line =~ s/\/\d+\-\d+\s*$//;	# >CP1C2Cuu3C4En5Na6Ar7con/1-762

			$line =~ s/\s\d+\.\.\d+$//;# pynast output puts these on: >uniques_1964 1..91


			if ($entry_counter >= 1)
				{

				$current_sequence =~ s/\s//g;$current_sequence =~ s/\t//g;
				if ($entry_counter == 1){$first_seq_length = length($current_sequence)};
				$current_seq_length=length($current_sequence);
				unless($first_seq_length == $current_seq_length)
					{print "\nWARNING ... unaligned .. first_seq_length:$first_seq_length, whereas $current_id has seq_length:$current_seq_length \n"};
#print "\ncurrent_id:$current_id current_seq_length:$current_seq_length current_sequence:$current_sequence\n";die;
				$seqs[$entry_counter] = $current_sequence;
				$alpha{$current_id} = $current_sequence;

	#			die"";

				};
			if ($line=~/^\s{0,2}>\s{0,2}(\S.*\S)\s*$/ )
				{
				$entry_counter ++;# print "ec:$entry_counter line:$line cid:$current_id\n";
				$current_id = $1; if($replace_id_spaces==1){ $current_id=~ s/\s/_/g};
				}else{print "error message! quitting";die};

			$ids[$entry_counter] = $current_id;
			$ids_r{$current_id} = 1;
			$current_sequence = "";
			}else{
			$current_sequence = $current_sequence . $line;
			if ($line=~/>/){print "WARNING: unreadable id line. quitting.$line\n";die};
			};
		}; #foreach $line(@all_lines)

	
	$current_sequence =~ s/\s//g;$current_sequence =~ s/\t//g;# bugfix 06MAR2015
	$ids[$entry_counter] = $current_id;
	$seqs[$entry_counter] = $current_sequence;
	$ids_r{$current_id} = 1;

print "current_id:$current_id\n";

#	close (FASTA_IN);

#if($screen_print == 1)
#	{
	print "$entry_counter entries read into memory. assuming all are length $current_seq_length\n";
#	};




};#sub read_fasta



####################################################################################
#
#
#
#
#
####################################################################################


		#### NO LONGER USED. SEE READ_MACCLADE INSTEAD
sub read_nexus
	{

if($screen_print == 1)
	{
	print "subroutine to read nexus format.\nreading input file.\n";
	};

	open(NEXUS_IN, $input) || die "Cant open input:$input.\n";
	

	my $in_matrix = 0;

	my $file_as_string = ""; my @all_lines = ();
	while($line= <NEXUS_IN>){$file_as_string .= $line};
	close(NEXUS_IN);
	@all_lines = split /\012\015?|\015\012?/, $file_as_string;

	foreach $line(@all_lines)
		{

#	while($line = <NEXUS_IN>)
#		{ #cb1
		
		$line =~ s/\n//;
		$line =~ s/\r//;

		if ($line =~ /\[Name\:\s+\S+\s+\S+\s*\S*\s*\S*\s*\S*\s*\S*\s*\S*\s*\S*\s*\S*\s+Len\:/i){print "warning: you have spaces in your ids, this isnt going to work. quitting\n";die}
		if ($line =~ /^\s*matrix\s*$/i){$in_matrix = 1};
		if ($line =~ /\;/ && $in_matrix == 1){$in_matrix = 0};
	#	print "$in_matrix\n";

		if ($in_matrix == 1)
			{
		#	 print "$line\n";
			if ($line =~ /^\s*(\S+)\s+(\S.+)$/)
				{
				my $current_id = $1;
				my $current_seq = $2;
				$current_seq =~ s/ //g;
				$seq_hash{$current_id} = $seq_hash{$current_id} . $current_seq;
				};

			};

		};
# close(NEXUS_IN);
	
my @key_list = keys %seq_hash;	print $#key_list , ", index of final key\n";

for $i(0 .. $#key_list)	{
my $current_key = $key_list[$i];$ids[$i+1] = $current_key;$seqs[$i+1] = $seq_hash{$current_key};
$entry_counter = $i+1;$current_seq_length = length($seqs[$i])
			};

	};




####################################################################################
#
#
#
#
#
####################################################################################

sub read_macclade
	{
	my $interleaved_counter = 0;
if($screen_print == 1)
	{
	print "subroutine to read macclade format.\nreading input file.\n";
	};

	open(MAC_IN, $input) || die "Cant open input:$input.\n";
	

	my $in_matrix = 0;

	my $file_as_string = ""; my @all_lines = ();
	while($line= <MAC_IN>){$file_as_string .= $line};


	$file_as_string =~ s/\[[^\]]+\]//g; # remove all comments

# open(OUT14, ">test_out") || die "cant open out";
# print OUT14 $file_as_string;
# close(OUT14);
# die;

	close(MAC_IN);
	@all_lines = split /\012\015?|\015\012?/, $file_as_string;

	foreach $line(@all_lines)
		{

# 'Lema_lichenis_BMNH*676188'                   TGCATGTCTCAGTA CAAGCC------ ----------   [2032]
	
		$line =~ s/\n//;
		$line =~ s/\r//;

		if ($line =~ /^\s*matrix\s*$/i){$in_matrix = 1};
		if ($line =~ /\;/ && $in_matrix == 1){$in_matrix = 0};
#		print "$in_matrix\n";

		if ($in_matrix == 1)
			{
#			 print "$line\n";

			if($line =~ /^(\'\s*\S.+\s.+\S\')\s+(\S+.+\S)\s*$/) # ebi readseq even puts spaces in sequence (why god?)
				{	# tax id includes spaces (wtF!)
				my $current_id = $1;
				my $current_seq = $2; $current_seq =~ s/\s//g;
				if(exists($seq_hash{$current_id})){$interleaved_counter++};
				$seq_hash{$current_id} .= $current_seq;
				
				}else{
				if ($line =~ /^\s*(\S+)\s+(\S+.+\S)\s*$/)
					{
					my $current_id = $1;
					my $current_seq = $2; $current_seq =~ s/\s//g;
					if(exists($seq_hash{$current_id})){$interleaved_counter++};
					# print "$current_id cs: $current_seq\n";
					$seq_hash{$current_id} .= $current_seq;
				
					};

				}

			};

		};

	
	my @key_list = keys %seq_hash;
	for $i(0 .. $#key_list)
		{
		my $current_key = $key_list[$i];$ids[$i+1] = $current_key;$seqs[$i+1] = $seq_hash{$current_key};
		$entry_counter = $i;$current_seq_length = length($seqs[$i]);
		};
	if($interleaved_counter>=1 && $interleaved_counter < $#key_list){$print_warning_message=1};
	};




####################################################################################
#
#
#
#
#
####################################################################################


sub read_clustal
	{

if($screen_print == 1)
	{
	print "subroutine to read clustal format.\nassumes no spaces or stars in matrix or id's\nreading input file.\n";
	};
	open(CLUSTAL_IN, $input) || die "Cant open input:$input.\n";
	

	my $in_matrix = 1;$i=0;

	my $file_as_string = ""; my @all_lines = ();
	while($line= <CLUSTAL_IN>){$file_as_string .= $line};
	close(CLUSTAL_IN);
	@all_lines = split /\012\015?|\015\012?/, $file_as_string;

	foreach $line(@all_lines)
		{

#	while($line = <CLUSTAL_IN>)
#		{ #cb1
		$line =~ s/\n//;
		$line =~ s/\r//;
		if ($line =~ /CLUSTAL|\*|multiple sequence alignment|MUSCLE/i){$in_matrix=0}else{$in_matrix=1}; 

			if ($line =~ /^(\S{2,80})\s{2,80}(\S+)$/ && $in_matrix == 1)
				{

# from tcoffe output:
#NCBI__Gastrancistrus_sp__D0730                                               CACGAGACCGATAGCGAACAAGTACCGTGAGGGAAAGTTGAAAAGAACTT


				# print "$line\n";

				my $current_id = $1;
				my $current_seq = $2;
				$current_seq =~ s/ //g;

				if(length($seq_hash{$current_id})>=1){}else{
					$i++;
					$ids[$i] = $current_id;
					};

				$seq_hash{$current_id} = $seq_hash{$current_id} . $current_seq;

				};



		};

# close(CLUSTAL_IN);

my @key_list = keys %seq_hash;

for $ii(1 .. $i)
	{
	my $current_key = $ids[$ii] ; $seqs[$ii] = $seq_hash{$current_key};
	$entry_counter = $ii;$current_seq_length = length($seqs[$ii])
	};

if($screen_print == 1)
	{
print "$entry_counter keys in seq hash\n";
	};
	};









####################################################################################
#
#
#
#
#
####################################################################################



sub read_phylip
	{

if($screen_print == 1)
	{
	print "subroutine to read phylip format.\nreading input file.\n";
	};

	open(PHYLIP_IN, $input) || die "Cant open input:$input.\n";
	

	my $in_matrix = 0;

	my $file_as_string = ""; my @all_lines = ();
	while($line= <PHYLIP_IN>){$file_as_string .= $line};
	close(PHYLIP_IN);

	unless($file_as_string =~ /^\s*\d+/){die "\n\nerror. you specified read phylip. but this doesnt look like phylip format.\n\n"}
	@all_lines = split /\012\015?|\015\012?/, $file_as_string;

	foreach $line(@all_lines)
		{
#	while($line = <PHYLIP_IN>)
#		{ #cb1
		$line =~ s/\n//;
		$line =~ s/\r//;
		if($line =~ /^\s*\d+\s+\d+\s*$/){$in_matrix = 0}else{$in_matrix = 1};

					# print "$line\n";
		if($in_matrix == 1)
			{	
			if ($line =~ /^\s*(\S+)\s+(.+)$/)
				{
				my $current_id = $1;#print "current_id:$current_id\n";
				my $current_seq = $2;
				$current_seq =~ s/ //g;
				$seq_hash{$current_id} = $seq_hash{$current_id} . $current_seq;
				push @key_list , $current_id;
				};
			};
		

		};
# close(PHYLIP_IN);

#my @key_list = keys %seq_hash;

for $i(0 .. $#key_list)
{my $current_key = $key_list[$i];$ids[$i+1] = $current_key;$seqs[$i+1] = $seq_hash{$current_key};$entry_counter = $i;$current_seq_length = length($seqs[$i])};
	};







####################################################################################
#
#
#
#
#
####################################################################################




sub write_nexus
	{

if($screen_print == 1)
	{
	print "subroutine to write nexus\nwriting to:$output\n";
	};

	open(NEXUS_OUT ,">$output") || die "cant open output file:$output\n";
	
	print NEXUS_OUT "\#NEXUS\n\nBEGIN DATA\;\nDIMENSIONS NTAX=$entry_counter NCHAR=$current_seq_length\;\n";
	print NEXUS_OUT "FORMAT DATATYPE=DNA INTERLEAVE GAP=- MISSING=?\;\n";


	if ($output_format =~ /nexus_concise/)
		{}else{
		for $i(1 .. $entry_counter)
			{
			$current_name = $ids[$i];
			print NEXUS_OUT "\[Name: ";

			if($current_name =~ /\S\s\S|\S-\S|\S\+\S|\S\(\S/) # 
				{
				print NEXUS_OUT "\'$current_name\'";
				}else{
				print NEXUS_OUT "$current_name";
				}

			print NEXUS_OUT "         Len:  $current_seq_length   Check:    0]\n";
			};
		};

	print NEXUS_OUT "\nMATRIX\n";

	for $i(1 .. $entry_counter)
		{
		$current_name = $ids[$i];
		$current_seq = $seqs[$i];
		if(length($current_name)<=1 || length($current_seq)<=1){print "warning: zero length of current entry. quitting\n";die}

		if($current_name =~ /\S\s\S|\S-\S|\S\+\S|\S\(\S/) # |\S-\S
			{
			print NEXUS_OUT "\'$current_name\'  $current_seq\n";
			}else{
			print NEXUS_OUT "$current_name  $current_seq\n";
			}
		};

	print NEXUS_OUT "\n\;\nEND\;\n";
	close(NEXUS_OUT);


	};



####################################################################################
#
#
#
#
#
####################################################################################



sub write_fasta
	{
my $count_printed_taxa=0;
if($screen_print == 1)
	{
	print "subroutine to write fasta\nwriting to:$output\n";
	};

	open(FASTA_OUT ,">$output") || die "cant open output file:$output\n";
	#open(FASTA_OUT_R ,">$output_r") || die "cant open output file:$output_r\n";
	#open(FASTA_OUT_A ,">$output_a") || die "cant open output file:$output_r\n";


	my @key_list_r = keys %ids_r;
#	fisher_yates_shuffle( \@key_list_r);
	@ids_a = sort @ids;

	for $i(1 .. ($entry_counter+1))# 04DEC2015: bugfix, last entry was previously ommited.
		{
		$current_name = $ids[$i]; # print "$i $current_name\n";

		my $current_name_r = $key_list_r[$i-1];
		my $current_name_a = $ids_a[$i];
		$current_seq = $seqs[$i];
		$current_seq_a = $alpha{$current_name_a};
		if(length($current_name)<=1 || length($current_seq)<=1){print "warning: zero length of current entry. quitting\n";die}
		if ($output_format =~ /fasta_degap/)
			{
			$current_seq =~ s/\-//g;
			$current_seq =~ s/\?//g;
			}
		if ($output_format =~ /fasta_missing_all/)
			{
			$current_seq =~ s/[\-\.]/$missing_data_character/g;
			}
		if ($output_format =~ /fasta_missing_terminal/)
			{
			if($current_seq =~ /^(\-+)([^\-].*)$/)
				{
				$the_gaps = $1;
				$the_rest = $2;
				$the_gaps =~ s/\-/$missing_data_character/g;
				$current_seq =~ s/^(\-+)([^\-].*)$/$the_gaps$the_rest/;
				};
			if($current_seq =~ /^(.*[^\-])(\-+)$/)
				{
				$the_gaps = $2;
				$the_rest = $1;
				$the_gaps =~ s/\-/$missing_data_character/g;
				$current_seq =~ s/^(.*[^\-])(\-+)$/$the_rest$the_gaps/;
				};

			}

		if($parse_taxa == 1)
			{
			foreach $taxa_parsee(@taxa_parse_list)
				{
				if($current_name =~ /$taxa_parsee/)
					{
					print FASTA_OUT ">$current_name\n$current_seq\n";
					print FASTA_OUT_R ">$current_name_r\n$current_seq\n";
					print FASTA_OUT_A ">$current_name_a\n$current_seq_a\n";
					$count_printed_taxa++;
					}
				}
			}else{
			if($current_seq =~ /^[\-N]+$/)
				{
				print "WARENING, removing entry with no sequence data ($current_name)\n";
				}else{
			print FASTA_OUT ">$current_name\n$current_seq\n";
		#	print FASTA_OUT_R ">$current_name_r\n$current_seq\n";
		#	print FASTA_OUT_A ">$current_name_a\n$current_seq_a\n";
			$count_printed_taxa++;

				}
			}
		};
	close(FASTA_OUT);
#	close(FASTA_OUT_R);
#	close(FASTA_OUT_A);

print "number printed to fastafile:$count_printed_taxa\n";

	};


sub write_phylip
	{

if($screen_print == 1)
	{
	print "subroutine to write phylip\nwriting to:$output\n";
	};
	my $should_i_print_warning = 0;
	open(PHYLIP_OUT ,">$output") || die "cant open output file:$output\n";
	
	print PHYLIP_OUT "$entry_counter $current_seq_length\n";

	for $i(1 .. $entry_counter)
		{
		$current_name = $ids[$i];
		$current_seq = $seqs[$i];

		if($remove_bad_bases == 1)
			{
			while($current_seq =~ s/J/?/){$bad_bases_removed++};
			};


		if(length($current_name)<=1 || length($current_seq)<=1){print "warning: zero length of current entry. quitting\n";die}
		$name_length = length($current_name);
	#	 print "current_name:$current_name name_length:$name_length\n";
		if($name_length == 1){print PHYLIP_OUT "$current_name         $current_seq\n"};
		if($name_length == 2){print PHYLIP_OUT "$current_name        $current_seq\n"};
		if($name_length == 3){print PHYLIP_OUT "$current_name       $current_seq\n"};
		if($name_length == 4){print PHYLIP_OUT "$current_name      $current_seq\n"};
		if($name_length == 5){print PHYLIP_OUT "$current_name     $current_seq\n"};
		if($name_length == 6){print PHYLIP_OUT "$current_name    $current_seq\n"};
		if($name_length == 7){print PHYLIP_OUT "$current_name   $current_seq\n"};
		if($name_length == 8){print PHYLIP_OUT "$current_name  $current_seq\n"};
		if($name_length >= 9){print PHYLIP_OUT "$current_name  $current_seq\n";$should_i_print_warning=1};

		if(exists($names_printed_to_phylip{$current_name}))
			{
			print "warning! you have duplicated sequence IDs ($current_name)\n";
			};
		$names_printed_to_phylip{$current_name} = 1;

		};
	close(PHYLIP_OUT);
	
	if($should_i_print_warning==1)	{
		print  "WARNING: phylip normally has 8 char limit on id. sequences have been converted nontheless\n";
		print  "as some programs (eg raxml ) accept relaxed phylip, without limitation on id length\n";
					};

	print "count entries printed to phylip outfile:$entry_counter
";
	};





####################################################################################
#
#
#
#
#
####################################################################################




sub write_paup
	{

if($screen_print == 1)
	{
	print "subroutine to write nexus\nwriting to:$output\n";
	};

	open(NEXUS_OUT ,">$output") || die "cant open output file:$output\n";
	
	print NEXUS_OUT "\#NEXUS\n\nBEGIN DATA\;\nDIMENSIONS NTAX=$entry_counter NCHAR=$current_seq_length\;\n";
	print NEXUS_OUT "FORMAT DATATYPE=DNA INTERLEAVE MISSING=-\;\n";


	if ($output_format =~ /nexus_concise/)
		{}else{
		for $i(1 .. $entry_counter)
			{
			$current_name = $ids[$i];
			print NEXUS_OUT "\[Name: $current_name         Len:  $current_seq_length   Check:    0]\n";
			};
		};

	print NEXUS_OUT "\nMATRIX\n";

	for $i(1 .. $entry_counter)
		{
		$current_name = $ids[$i];
		$current_seq = $seqs[$i];
		if(length($current_name)<=1 || length($current_seq)<=1){print "warning: zero length of current entry. quitting\n";die}

		print NEXUS_OUT "$current_name  $current_seq\n";
		};

	print NEXUS_OUT "\n\;\nEND\;\n";

if ($output_format =~ /paup_bootstrap/i)
	{

	print NEXUS_OUT "
set autoclose = yes\;
set storetreewts=yes\;
set maxtrees=10000\;
log file=paup_logfile\;\n";

if(length($outgroup_name)>=2)
	{print NEXUS_OUT "OUTGROUP $outgroup_name\;\n";
	};

print NEXUS_OUT "bootstrap nreps=10 search=heuristic treefile=bootstrap_treefile1 replace=no brlen=yes/ start=stepwise addseq=random nreps=4 savereps=no randomize=addseq hold=1 swap=tbr multrees=yes\;
savetrees from=1 to=1 file=savetrees_file format=altnexus brlens=yes savebootp=nodelabels MaxDecimals=0 \;
contree all /majrule=yes strict=no le50=yes usetreewts=no showtree=yes treefile=contree_file replace=yes grpfreq=yes\;
log stop\;
\[to run type: paup input_file.nex -n\]\n";

	};

if ($output_format =~ /paup_dist/i)
	{
	print NEXUS_OUT "\nBEGIN PAUP\;\n set autoclose = yes\;\n";
	if(length($outgroup_name)>=2)
		{print NEXUS_OUT "OUTGROUP $outgroup_name\;\n";
		}
	print NEXUS_OUT "DSET DISTANCE=GTR RATES=GAMMA NEGBRLEN=PROHIBIT\;\n";
	print NEXUS_OUT "NJ BRLENS=YES BREAKTIES=RANDOM SHOWTREE=YES TREEFILE=distance_treefile\;\n";
	print NEXUS_OUT "\[to run type: paup input_file.nex -n\]\n";
	};

if ($output_format =~ /paup_pars/i)
	{
	print NEXUS_OUT "\nBEGIN PAUP\;\n set autoclose = yes\;\n";
	if(length($outgroup_name)>=2)
		{print NEXUS_OUT "OUTGROUP $outgroup_name\;\n";
		};
	print NEXUS_OUT "hsearch multrees=no addseq=random nreps=100\;\n";
	print NEXUS_OUT "showtrees 1 /userbrlens=no\;\n";
	print NEXUS_OUT "PSCORES 1 /ci=yes ri=yes tl=yes\;\n";
	print NEXUS_OUT "\[to run type: paup input_file.nex -n\]\n";
	};


	print NEXUS_OUT "\nQUIT\;\nENDBLOCK\;\n";
	close(NEXUS_OUT);


	};



####################################################################################
#
#
#
#
#
####################################################################################




