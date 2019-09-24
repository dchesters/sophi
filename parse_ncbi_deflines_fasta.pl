
###################################################################################################################################
#
#
#
#	parse_ncbi_deflines_fasta.pl
#
#    	Copyright (C) 2016-2017  Douglas Chesters
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
#	parses fasta files in the format present when directly downloaded from genbank.
#	(at least, the proteome data i was using)
#	
#	
#	
#	
#	
#	
#	change log
#	08 Jan 2015: adapted for deflines of transcriptome entries (in the WGS division).
#	15 Jan 2015 (1.02): bufix reverting to proteomes
#	03 Oct 2015: reads lines with some NP_ number not GI
#	30 Oct 2015: if no NP_ nor GI, just makes up some meaningless number to use instead
#	06 Apr 2016: setting given in command line,
#		so more user friendly when in a pipeline (user doesnt need to open script and change setting)
#	25 Apr 2016: parse errors are printed to a log file, since i run this in a big loop, 
#		otherwise would be easy to miss this happened
#	12 Jul 2017: print parsed taxon name to the log file
#	17 Aug 2019: resurected! ncbi long since stopped using gi numbers, minor change for using accession.
#	
#	
#	
#	
#	
#	
#	
#################################################################################################





$arg_string  = join ' ', @ARGV;

#####################################
read_command_arguments($arg_string);#
#####################################



$verbose = 0;

if($id_format == 1)
	{print "\nyou have chosen \$id_format == 1: printing fasta IDs >[taxon]_[GI]\n"};
if($id_format == 2)
	{print "\nyou have chosen \$id_format == 2: printing fasta IDs >[GI]\n"};





open(IN, $in) || die "\nerror cant open nin $in\n";
open(OUT, ">$in.b") || die "\nerror cant open out\n";

if($arg_string =~ /-output_key/)
	{open(KEY, ">$in.key") || die "\nerror cant open out\n"};

my $parsed_taxon = "NA";

while ( my $line = <IN>)
	{

	my $original_line = $line;
	# refseq proteome defline looks like this:
	#>gi|328696447|ref|XP_003240026.1| PREDICTED: uncharacterized protein LOC100569617 [Acyrthosiphon pisum]
	# transscriptome deflines look like this:
	#gi|712074138|gb|GASQ01000001.1| TSA: Tetrix subulata s1_L_1_0 transcribed RNA sequence

	# C. elegans proteome fasta file IDs look like this, no GI, parse that NP number instead:
	#>NP_001020953.1 Phosphatidylinositol 3-kinase [Caenorhabditis elegans]
	#  turned out a little awkward to parse as with others, so quick solution:
	$line =~ s/^\>NP_(\d+)\.\d\s/">gi|" . $1 . "|string"/e;


	# if no GI number , make one:
#	unless($line =~ />gi\|/)
#		{
#		$line =~ s/^>/>gi|9999$meaningless_number|/;
#		$meaningless_number++;
#		}

	if($line =~ />([^\.\s]{6,10})(.+)/)
		{
		my $gi = $1;my $rest = $2;
#		}elsif($line =~ /^\>(NP_\d+\.\d)\s(.+)/)
		#print $line;die;
		$count_entries++;
		 my $tax = "NA";
		$rest =~ s/\s([A-Z][a-z]+)\ssp\.\s/ $1 sp /;# get rid of dot after sp in the unidentified species
		$rest =~ s/(\[[A-Z][a-z]+)\ssp\.(\])/$1 sp$2/;# get rid of dot after sp in the unidentified species

		$rest =~ s/ ([A-Z][a-z]+)\snr\. ([a-z]+)\s/ $1 $2 /; #  Contacyphon nr. frater C


		#print "rest:$rest\n";
		if($rest =~ /\[([A-Z][a-z]+\s[a-z]+)\]/)
			{
			$tax = $1
			}elsif($rest =~ /\[([A-Z][a-z]+\s[a-z]+)\s[a-z]+\]/)# subsp: [Apis mellifera ligustica]
				{
				$tax = $1
				}elsif($rest =~ /\s([A-Z][a-z]+\s[a-z]+)\s/)
					{
					$tax = $1
					}else{
					print "line:($line)\ngi:($gi) rest:($rest)\n";
					die "\nerror 138. cant parse binomial from this string:($rest)\n";
					};



		unless($tax =~ s/^([A-Z][a-z]+\s[a-z]+).*/$1/)
			{die "\nerror, unexpected taxon name:($tax)\n"}
		$tax =~ s/\s+/_/; $parsed_taxon = $tax;

		if($verbose == 1)
			{
			print "gi:$gi tax:$tax\n";
			};			

		my $new_ID = "NA";
		if($id_format == 1){	$new_ID = "$tax" . "_$gi"};
		if($id_format == 2){	$new_ID = "$gi"};
		if(exists($all_new_IDs{$new_ID}))
			{print "\nwarning. after reformatting ids, this one is duplicated:$new_ID\n"}; 
		$all_new_IDs{$new_ID} =1;
		print OUT ">$new_ID\n";
		print KEY "$new_ID\t$original_line";

		if($count_entries =~ /0000$/)
			{print "fas $count_entries, gi:$gi tax:$tax, newID:$new_ID\n"}


		}else{
		if($line =~ />/)
			{
			print "line:$line";$date = localtime time;
			open(ERROR_LOG, ">>parse_ncbi_deflines_ERROR_LOG");
			print ERROR_LOG "parse error for file $in, parsing incomplete after $count_entries entries at time $date\n";
			close ERROR_LOG;
			die "\n\nerror 14 . parse error, ID line doesnt match regex />gi\|([^\|]+)\|(.+)/\n$line\n";
			};

		print OUT "$line";
		};

	};
close IN;
close OUT;

my @ids = keys %all_new_IDs;
print scalar @ids, " new ids
count_entries:$count_entries

script has modified fasta IDs, if you want to know the original ids, use option -output_key
\n\n";



open(LOG , ">>Parse_NCBI_deflines_LOG");
print LOG "$count_entries\t$in\t$parsed_taxon\n";
close LOG;

close KEY;





#########################################################################################################

sub read_command_arguments
{

my $arguments = shift;
$arguments =~ s/\n//;$arguments =~ s/\r//;


# $in = $ARGV[0];


# $id_format = 1;# 1= >taxon_GI; 2= >GI

if($arguments =~ /-in\s+(\S+)/)
	{
	$in = $1;
	open(INTEST , $in) || die "\ninput not found in current directlry:$in\n";
	close INTEST;
	}else{
	print "\nerror reading command arguments  (-in)\n";die""
	}

if($arguments =~ /-format\s+([12])/)
	{
	$id_format = $1;
	}else{
	print "\nerror reading command arguments  (-format)\n";die""
		
	};

};

#########################################################################################################








