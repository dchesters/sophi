
my $db 				= $ARGV[0];
my $uclust_outfile 		= $ARGV[1];
my $sequence_identity_cutoff 	= $ARGV[2];
my $trim_retrieved_sequences 	= $ARGV[3];	# 0= default, no trim; 1=trim according to positions of longest hit; 
						# 2=trim according to left-most and right-most positions from all hits to given database sequence
my $sequence_length_limit 	= $ARGV[4];
my $parse_which		 	= $ARGV[5];	# column 1 or 2, query or hit. default is hit (2)
my $blastdbcmd 			= $ARGV[6];	#"blastdbcmd-2.2.28+-64bit";# 2.2.25 not installed on thinkpad



###################################################################################################################################
#
#
#
# 	parse_hits.pl, Perl script for parsing output of blast or uclust search
#
#    	Copyright (C) 2013-2014 Douglas Chesters
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
# 	to run:
# 	perl parse_hits.pl database.fas blast_output 25
# 	the third argument (the number) is percent identity cutoff. hits with identity below this value will be ignored.
#
# 	input is blast or uclust search results in tabular format .... command in blast to do this is:
# 	makeblastdb -in db.fas -dbtype nucl -parse_seqids
# 	blastn -task blastn -db db.fas -out blastoutTEST -dust no -strand both -evalue 1e-6 -query queries.fas -num_threads 1 -num_descriptions 100000 -num_alignments 100000 -outfmt '6 qseqid sseqid evalue pident length sstart send qframe sframe'
#
# 	and the uclust command for older version (4.2.66) is:
# 	usearch --query $query_file --db db.fas --maxlen 50000 --userout uclust_out --userfields query+target+evalue+id+cols+tlo+thi+ts+strand --evalue 1e-6 --nousort
#
# 	for newer version (6.0.3) should be:
# 	usearch -search_global queries.fas -db db.fas -strand both -id 0.25 -userout AllplantSeqsEdited.AAA -userfields query+target+evalue+id+alnlen+tlo+thi+ts+qstrand
#
#	
#	
# 	change log:
# 	16 sep 2013, length limit applied to output seqs in addiiton to aligned segment.
# 	oct 2013: GPL added
#	9 oct 2013: sequence length limit given by user in command line, as last argument
#	7 mar 2014: option to parse column 1 or 2, usually query followed by hit
#	12 jun 2014: bugfix, detects presence of blast searchable db where it is split across multiple files (if (-e "$db.00.nhr"))
# 	26 Aug 2014: blastdbcmd given in command line
#	04 Jan 2015: autodetect protein blast database also
#		will extract via gi numbers, if database is using standard ncbi deflines
#	
#	
#	
#	
#	
#
#
# 	********** USER OPTIONS BELOW **********
#


my $verbose			= 0;	# if = 1, for debugging, prints loads of details to screen, but will slow it down a lot.
my $print_warnings		= 0;	# if = 1, warns of sequences in reverse orientation but not minus strand, or other way round.

# 	if the database specified by the user is blast compatible (via makeblastb), the script will use blastdbcmd to get sequences. 
# 	if it is not, the script will try to read the whole db into memory.
# 	under what name did you install blastdbcmd? write this name into the following variable:


my $blastdb_present 		= 0; 	# default == 1, retrieves hits using blast database 
					# plus the software blastdbcmd
					# if for whatever reason these are not available,
					# select option 0 here, 
					# then this script itself reads the database for retrieval of hits. 
					# note sequence trimming is not supported for == 0

#
# 	********** END OF USER OPTIONS **********
#
#
#####################################################################################################################################################



$extract_using_gi 	= 0;# 0=for my traditional IDs. 1=genbank deflines containing GI's

$dbtype 		= "nucl"; # "nucl" or "prot"



##################
check_settings();#
##################


my %lengths = ();
my $current_hit_strand =0;
my $current_hit_reversed =0;
%db_stored  =();

unless($sequence_length_limit =~ /\d/){$sequence_length_limit = 50}

open(IN, $uclust_outfile) || die "\n\nerror cant open file:$uclust_outfile\n";

print "\nparse_uclust_searchresults.pl file:$uclust_outfile db:$db sequence_length_limit:$sequence_length_limit\n\n";

my $count_forward_hits =0;my $count_reverse_hits =0;
my $count_plus_strand =0;my $count_minus_strand =0;

my $linecount=0;my $count_short_hit_discarded = 0;

while ( my $line =<IN>)
	{
	$line =~ s/\n//;$line =~ s/\r//;
#gi|642930308|ref|XP_008196340.1|	1	835	835	0.0	100.00	1

	if($line =~ /^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
		{
		my $query= $1; my $subject = $2;my $pident = $4;
		my 		$length=$5; my 	$tlo=$6; my 	$thi=$7;my 		$strand = $9;
		# uclust query+ target+evalue+id+	cols+		tlo+		thi+		ts+	strand
		
		# blast: qseqid sseqid evalue pident 	length 		sstart 		send 		slen 	sframe



		my $target;
		if($parse_which == 1){$target = $query}else{$target = $subject}

		$linecount++;
		if($verbose == 1)
			{
			print "$uclust_outfile line $linecount\n"
			}else{
			if($linecount =~ /00000$/){print "$uclust_outfile line $linecount\n"}
			}

		if($length <= $sequence_length_limit)
			{
			$count_short_hit_discarded++;
			}else{

		my $left = "";my $right = "";

		#unless($linecount == 1)# why was this here, header previously?
		#	{

			my $reversed = 0;
			if($tlo <= $thi )
				{
				$left = $tlo; $right = $thi;
				#$reversed{$target} =0;
				$current_hit_reversed =0;
				$count_forward_hits++;if($verbose == 1){print "not reversed\n"}
				}else{
				$left = $thi; $right = $tlo;
				#$reversed{$target} =1;
				$current_hit_reversed =1;
				$count_reverse_hits++;if($verbose == 1){print "is reversed\n"}
				}

			if($strand =~ /[\+1]/)
				{
				$current_hit_strand="plus";
				#$strand{$target}="plus";
				$count_plus_strand++;
				}else{
				$current_hit_strand="minus";
				#$strand{$target}="minus";
				$count_minus_strand++;
				}

			if($print_warnings == 1)
				{
				if($current_hit_reversed == 1 && $current_hit_strand eq "plus")
					{print "Warning1: hit ($target) appears reversed (tlo:$tlo > thi:$thi) although strand ($strand) is plus \n$line\n\n"}
				if($current_hit_reversed == 0 && $current_hit_strand eq "minus")
					{print "Warning2: hit ($target) appears not reversed (tlo:$tlo < thi:$thi) although strand ($strand) is minus\n"}
				}

			#};

		unless($left =~ /\d/ &&  $right=~ /\d/ ){die "\nerror\n"}

		#if($line =~ /642930308/){print "q$query s$subject p$pident \n"};

		if(exists($start{$target}))
			{
			# instead of taking left-most and right-most, which seems v liberal, take longest hit
			#if($line =~ /642930308/){print "p$pident, s$sequence_identity_cutoff, l$left, r$right, $length\n"};
			if($pident >= $sequence_identity_cutoff)
				{
				if($verbose == 1){print "$target has been observed already. identity ($pident) >= cutoff ($sequence_identity_cutoff) ... comparing to previous.\n"}

				##################################################
				store_positions($left, $right, $length, $target);#
				##################################################

				}else{

				if($verbose == 1){print "$target has been observed already. identity ($pident) < cutoff ($sequence_identity_cutoff) ... so ignoring.\n"}
				}
			}else{
			#if($line =~ /642930308/)
			#	{
			#	print "p2:$pident, s2:$sequence_identity_cutoff, l2:$left, r2:$right, $length\n";
			#	print "$line\n"
			#	};
			if($pident >= $sequence_identity_cutoff)
				{
				if($verbose == 1){print "$target has NOT been observed already, and identity ($pident) >= cutoff ($sequence_identity_cutoff), storing.\n"}
				#store_positions($left, $right, $target);
				$start{$target}=$left;$end{$target}=$right;$lengths{$target} = $length;
				$strand{$target}=$current_hit_strand; $reversed{$target} = $current_hit_reversed;

				}else{
				if($verbose == 1){print "$target has NOT been observed already, but identity ($pident) < cutoff ($sequence_identity_cutoff), NOT storing.\n"}
				}			
			}
			}


		
		}else{

		#print "test45:$line\n";
		}

	}

close(IN);




if($blastdb_present == 0)
	{read_db_into_memory();
}else{
open(OUTLOG, ">$uclust_outfile.retreivedINFO") || die "\nerror\n";
}

#die;


unless($linecount >= 2){die "\n\nerror. lines in search results file ($uclust_outfile) do not match regex\n"}

if($count_minus_strand == 0){print "\nnote: according to column 9 in the file:$uclust_outfile, there are no hits with minus strand. is this expected?\n"}

print "\nblast output line count ($linecount). short hits discarded ($count_short_hit_discarded). \nextracting seqs\n";

my @hits = keys %start;@hits = sort @hits;
my $tot = $#hits;
my $counted_done=0;
my $confirm_reversed =0;
my $confirm_minus_strand =0;


open(OUT, ">$uclust_outfile.retreived") || die "\nerror\n";
print OUTLOG "target start end reversed strand len\n";


foreach my $target(@hits)
{

my $len = $end{$target}-$start{$target};
#print "target:$target en:$end{$target} st:$start{$target}\n";
unless($end{$target} =~ /\d/){die "\nerror. no start or end positions parsed for entry ($target)\n"}
print OUTLOG "$target $start{$target} $end{$target} $reversed{$target} $strand{$target} $len\n";

my $trim_string = "";
unless($trim_retrieved_sequences == 0)
	{$trim_string = " -range $start{$target}-$end{$target}"};

#my $entry_retrieved =`blastdbcmd-2.2.25+-64bit -db inv.fas.parsed.rr -dbtype nucl -strand plus -line_length 60 -range 1-345 -entry ICP1Ch3C4C5C6Or7spec_GQ392362`;

# should be precisely nucl or prot

#-strand $strand{$target} ### strand option gets reverse complement, not just other strand

my $target2 = $target;

#print "target2 A:$target2\n";
if($extract_using_gi == 1){$target2 =~ s/gi\|(\d+)\|.+/$1/}#gi|642930308|ref|XP_008196340.1|
#print "target2 B:$target2\n";

my $sequence_extraction_string = "$blastdbcmd -db $db -dbtype $dbtype -line_length 100$trim_string -outfmt=\%s -entry $target2";

if($verbose == 1){
print "command to extract current seq is:$sequence_extraction_string\n";
}

#print "sequence_extraction_string:$sequence_extraction_string\n";

my $entry_retrieved = "";

if($blastdb_present == 1)
	{
	####################################################
	$entry_retrieved =`$sequence_extraction_string`;#
	####################################################
	}else{

	$entry_retrieved = $db_stored{$target};

	#print "\ntarget:$target\n\tseq:" , $db_stored{$target} ;

	}


#$entry_retrieved =~ s/>lcl\|/>/;
#$entry_retrieved =~ s/\:$start{$target}\-$end{$target}//;
#$entry_retrieved =~ s/\:\d+.\d+\s+No definition line found//;
#$entry_retrieved =~ s/\s*No definition line found//;
$entry_retrieved =~ s/\n//;$entry_retrieved =~ s/\r//;

if($reversed{$target}==1)
	{
	#print "entry1:$entry_retrieved\n";
	$entry_retrieved = reverse($entry_retrieved);
	$confirm_reversed++;
	#print "entry2:$entry_retrieved\n";

#	}

#if($strand{$target} eq "minus")
#	{
#	$confirm_minus_strand++;
#	#print "1 $target minus $entry_retrieved\n";
	$entry_retrieved =~ tr/ACGTYRacgtyr/TGCARYtgcary/;
#	#print "2 $target minus $entry_retrieved\n";

	}


if($entry_retrieved =~ /[actgATCG]/)
	{

	unless(length($entry_retrieved) <= $sequence_length_limit)########### sept 2013
		{
		if($verbose == 1){print "$counted_done of $tot. extracting sequence with $blastdbcmd, -entry $target -range $start{$target}-$end{$target} -strand $strand{$target} reversed:$reversed{$target}\n"};
		print OUT ">$target\n$entry_retrieved\n";$counted_done++;
		}

	}else{
	print "\n\nWARNING: no sequence retrieved for $target. if this is rare, can be ignored. \n";
	print "this problem often due to use of atypical characters in taxon labels.\n";
	print "another cause could be if you made some manual changes to the fasta database ($db), but did not remake the blastdb (makeblastdb).\n\n";
#	die;
	}

# old version of blast
#my $entry_retrieved =`fastacmd -L $current_start,$current_end -d endopterygota_fasta_coded -s $current_id -S $current_strand`;
	

if($counted_done =~ /000$/){print "$counted_done out of $tot\n"}



}

close(OUT);
close(OUTLOG);



print "\nin the blast/uclust output file:
forward_hits:($count_forward_hits) reversed_hits:($count_reverse_hits)
plus_strand:($count_plus_strand) minus_strand:($count_minus_strand)\n";

print "\n\n$linecount in $uclust_outfile. whereas $counted_done sequences have been printed to the outfile $uclust_outfile.retreived
remember, where multiple queries are used, there will be multiple hits to the same database entries.
\n
of the $counted_done sequences printed to the output file, $confirm_reversed were reversed and swapped strand, prior to printing. 
";


print "\nend of script.\n";
exit;


#############################################################################################################################################
#
#
#
#
#
#############################################################################################################################################



sub check_settings
{

unless($sequence_identity_cutoff =~ /\d/)
	{$sequence_identity_cutoff = 10; print "\nno sequence identity cutoff given in command. using default of 10 percent\n"}

unless ($db =~ /[a-z0-9]/i && $uclust_outfile =~ /[a-z0-9]/i)
	{die "\n\n\terror. to run this program, 
	specify both the name of the fasta database (from which the sequeneces will be obtained), and the blast / uclust output: 
	perl parse_uclust_searchresults.pl database_filename search_output_filename
	\n"}

unless($blastdbcmd =~ /blastdbcmd/)
	{
	die "\nerror, expecting blastdbcmd command as last argument (as of Aug 2014). quitting.\n";
	}



my $test2 =`$blastdbcmd -help`;
unless($test2 =~ /DESCRIPTION/)
	{print "\n\n\nerror. no response to the command $blastdbcmd, 
	check that what you have in the variable \$blastdbcmd (at the top of the parse_uclust_searchresults.pl script)
	 matches the name under which you installed that program.\n\n"}


if (-e "$db.nhr" || -e "$db.phr")
	{
	$blastdb_present=1; print "\nblast searchable db found\n"
	}else{
	if (-e "$db.00.nhr" || -e "$db.00.phr")
		{
		$blastdb_present=1; print "\nblast searchable db found\n"
		
		}else{
		$blastdb_present =0;print "\n\nwarning. $db is not blast searchable (i.e. via makeblastdb). it will be read into memory, trimming currently unsupported in such case.\n"
		}
	}


if($trim_retrieved_sequences == 0)
	{
	print "user set \$trim_retrieved_sequences to 0, not trimming hits\n"

	}elsif($trim_retrieved_sequences == 1)
	{
	print "user set \$trim_retrieved_sequences to 1, trimming according to positions of longest hit\n"

	}elsif($trim_retrieved_sequences == 2)
	{
	print "user set \$trim_retrieved_sequences to 2, trimming according to left-most and right-most positions from all hits to given database sequence\n"

	}else{
	print "user did not set \$trim_retrieved_sequences. using default, no trimming\n";$trim_retrieved_sequences=0;
	};


}




#############################################################################################################################################
#
#
#
#
#
#############################################################################################################################################




sub store_positions
{
my $left_position = $_[0];my $right_position = $_[1];my $the_length = $_[2];my $the_target = $_[3];


if($trim_retrieved_sequences =~ /[01]/)
	{

	if($the_length >= $lengths{$the_target})
		{# store new hit if its longer than previously stored hit
		$start{$the_target}=$left_position;$end{$the_target}=$right_position;$lengths{$the_target} = $the_length;
		$strand{$the_target}=$current_hit_strand; $reversed{$the_target} = $current_hit_reversed;
		if($verbose == 1){print "new length ($the_length) >= previous ($lengths{$the_target}) ... storing new\n"}
		}else{
		if($verbose == 1){print "new length ($the_length) < previous ($lengths{$the_target})  ... keeping old\n"}
		
		}

	}
if($trim_retrieved_sequences == 2)
	{


	if($left_position <= $start{$the_target})	
		{# store new position (left OR right) if its outside the range of previously stored
		$start{$the_target}=$left_position;
		$strand{$the_target}=$current_hit_strand; $reversed{$the_target} = $current_hit_reversed;
		if($verbose == 1){print "current left ($left_position) is more to left than stored ($start{$the_target}), adopting new\n"};
		}else{
		if($verbose == 1){print "current left ($left_position) is not more to left than stored ($start{$the_target}), keeping old\n"};

		};
	if($right_position >= $end{$the_target})
		{
		$end{$the_target}=$right_position;
		$strand{$the_target}=$current_hit_strand; $reversed{$the_target} = $current_hit_reversed;
		if($verbose == 1){print "current right ($right_position) more to right than stored ($end{$the_target}), adopting new\n"};

		}else{
		if($verbose == 1){print "current right ($right_position) not more to right than stored ($end{$the_target}), keeping old\n"};
		
		};

	}




}# sub store_positions





#############################################################################################################################################
#
#
#
#
#
#############################################################################################################################################


sub read_db_into_memory
{
print "sub read_db_into_memory\n";

my $db_as_string = "";
open(DB , "$db") || die "\nerror 432. cant open $db\n";
while (my $line = <DB>)	{$db_as_string .= $line};close(DB);
my @db_array = split />/, $db_as_string;
print scalar @db_array , " seqs in file\n";

for my $each_line(1 .. $#db_array)
	{
	my $line = $db_array[$each_line];
	if($line =~ /^(.+)/)
		{
		my $speciesid = $1;	#print "$speciesid\n";
		$line =~ s/^.+\n//;
		$line =~ s/\012\015?|\015\012?//g;$line =~ s/\n//g;$line =~ s/\r//g;
		$line =~ s/[\s\t]//g;

		$db_stored{$speciesid} = $line;

		}
	}




}



#############################################################################################################################################
#
#
#
#
#
#############################################################################################################################################


