











# 
# based on parse_ortholog_results.pl
# 
# -in [] -out_prefix [] -filter_duplicates -output_id_format
# 
# 
# 
# 
# queries:Buchnera protein IDs; hits:GInumber followed by arbitrary protein sequence ID
# need in output: each protein in a different file, fasta ids are GI numbers.
# 
# 
# 
# 
# 
# 
# what to filter????
# 	there are multiple queries, so more than one query would hit a subject, so filter hits to same thing.
# 	blast may report more than one sequence for a given hit (different HSP's)
# 	and more than one locus in a given transcriptome may be hit
# 	HSP's are presumably difficult to processm so just take the best one (longest)
# 	phylotreepruner will decide which locus is not a paralog, so dont filter different loci
# 	so we only want to filter so for each transcriptome we have one sequence per transcript
# 
# 	warning, current setup, must use phylotreepruner, there are multiple transcipts for each gene-species
# 
# 
# 
# 
# 
# change log
# 15 aug 2016: allowing multiple alignments per query, need to pick best one
# 18 may 2017: bugfix due to change in input format numbers for genes, to abreviated gene names
# 31 jul 2017: wont give error message for 1KITE's ortholog IDs 
# 
# 
# 
# 
# 
# 
# 










$arg_string  = join ' ', @ARGV;

#####################################
read_command_arguments($arg_string);#
#####################################




				# these are what are tested for duplication:
# $dont_print_duplicates = 1; # $query_protein_ID.binomial.transcript_number
				# i.e. you really dont want these.

$speciesnumber = 1;

$ortholog_taxon_limit = 47; # default 1, 


# $format = 3; 	# 1=species name ; 
		# 2= Gen_sp|code; arbitrary code for each species, for phylotreepruner
		# 3= Gen_sp|transcriptNumber; species_name transcript_number, for phylotreepruner (more suitable than above)

unless($string =~ /\w/){die "\nsyntax error. no string given \n"};

if($format == 1){print "\nyou selceted to print to outputs format 1. This is species name only\n"};
if($format == 2){print "\nyou selceted to print to outputs format 2. This is arbitrary code for each species, for phylotreepruner\n"};
if($format == 3){print "\nyou selceted to print to outputs format 3. This is species_name transcript_number, for phylotreepruner (more suitable than option 2)\n"};



# $annotation_info_file = "/home/douglas/usr_files/annotated_reference_sequences/insecta_hmmer3-2/insecta_hmmer3-2/annotation/insecta_hmmer3-2.annotation.ixosc.txt";
# open(FILE , "$annotation_info_file");
# while (my $line  =<FILE>)
# 	{
# 	if($line =~ /^(\d+)\|[^\|]+\|(.+)/)
# 		{my $code = $1; my $annot = $2;#print "code:$code annot:$annot\n";
# 		$store_annotations{$code} = $annot};
# 	}
# close FILE;



# first go through file, decide which lines to use...

open(IN, $in) || die "\nerror cant open in ($in)\n";
while (my $line = <IN>)
	{

	$line =~ s/\n//;$line =~ s/\r//;

	#qseqid sseqid evalue pident length sseq
	#646724193	Occasjapyx_japonicus_999973103	2e-49	48.04	179	DDCDAVVERILSDLDCAFRRDPLIDEFGIVPCYQRQNKSPVVHVEHKLGLEDWSVKPLYMYAYGKIMEWKRTR-RREEPDRLLQWTRAALLINPDITTFWNVRKTLASSGNLSVDAELHLTTLVLTRKPKSVETMSHRRWILGQ-SDSCSSEEDLLRTEMSVALASAGRYPNNYHAWSH
	# 

	# to be awkard, some have GI's, some do not
	#646724120	Ctenocephalides_felis_707807658	0.0	50.79	697	MAKLNKSFIWGTLFASLTWTVSLYLYWGLNSTAETQQSTSISYSPVPAKFIGKQRPYSSNDVL---AQNKNDVSRKALKWEKYAKSKYMLEERGMIKYDVL------KKPKMSEDMEKKLY-NPNKVSDKLLKELQPIVVLPEFDEPGMVHTIDEQILKEEGYKNHAFNVLVSKKIGLERNVPDTRHELCEKQVYDSNLPKASIVICFYNEDLTTLLRSVYSVFQKTPEHLIHEVILVDDYSDRCDELGEKLESDLRFMFEDSKGTLKPEWYKKVHILRMRNREGLIRARVLGARNATGEALVFLDSHIEVNKEWLQPLLQRIKQNETTVVMPVIDIINADTFNYTPSPLVRGGFNWGLHFKWDNLPKGTLSQDADFIKPIRSPTMAGGLFAINKDYFTHLGEYDMGMDIWGGENLEISFRIWMCGGSLELIPCSRVGHVFRKRRPYGSPGGTDSMLKNSLRVAHVWMDDYKDYFLKQHNQAKYVSYGNVQPRIKLRKELGCKSFKWYLENIYPELNLPGE-----KSTGKSDGEIKFQPWHSRKRDYVSNYQIRLQNTSYCLQSSERTGGQKEILKRGMHLTVQPCIRVKHQMWYETSKSELVLAQLLCLEASLGSSASGIGAAKSRPILGKCHEMGGDQEWKHKRENDTPIYNLASGTCLGVKRLPTTDSLAGAVVEMVLCSASGDEYSLRWDLV

# using core orthologs with 6 members for each gene, using protein id so same query ID 
#411975	Acerentomon_sp_9999172743	3e-34	35.29	255	FDSNKQMVAAGDAFPHKYCTKLDRGDYVIRLQIRHERRELLEGFTDCSLLLRQKLASAISLDVYSNMNNALSQGKKLSSMLLPHRHTAPLYIAA-PPQDKI--GKLIPPGSFFVATITLCKNEQHKKVDSYTIRLLLADPSKKASNLAKSDKKEAKNKVEEYDEASRDFKVTWLSKLDP--------------------NGVVGVKLYDELISCYPSHLPSLVARLQALDSDNRDLDNPSSKI--DKMDLTEEII
#411975	Acerentomon_sp_9999172743	6e-31	38.20	233	FDSNKQMVAAGDAFPHKYCTKLDRGDYVIRLQIRHERRELLEGFTDCSLLLRQKLASAISLDVYSNMNNALSQGKKLSSMLLPHRHTAPLYIAAPPQDKIGKLIPPGSFFVATITLCKNEQHKKVDSYTIRLLLADPSKKASNLAKSDKK--------------------EAKNKVEEYDEASRDFKVTWLSKLDPNGVVGVKLYDELISCYPSHLPSLVARLQALDSDNRDL
#411975	Acerentomon_sp_9999172743	3e-36	37.44	227	FDSNKQMVAAGDAFPHKYCTKLDRGDYVIRLQIRHERRELLEGFTDCSLLLRQKLASAISLDVYSNMNNALSQGKKLSSMLLPHRHTAPLYIAAPPQDKIGKL-IPPGSFFVATITLCKNEQHKKVDSYTIRLLLADPSKKASNLAKSDKKEAKNKVEEYDEASRDFKVTWLSKLDPNGVVGVKLYDELISCYPSHLPSLVARLQALDSDNRDLDNPSSKIDKMDLT
#411976	Acerentomon_sp_9999213752	0.0	66.59	910	KTLAIRREEQSIWERRAPLSPTHVRKLVKAGVRVIVQPSNRRAYPMQTYANAGAIVKEDISEAPVVFGVKQVPVDALQPNKTYCFFSHTIKAQEGNMPLLDAILKKNIRLIDYEKMVVQSGQRIVAFGKYAGVAGMINILHGLGLRLLALGHHTPFMHIGPAHNYRNSMMARQAVRDAGYQIALGLLPRSIGPLTIVFTGTGNVSQGAQEIFQELPHEYIPVEMLQKVAEHGATNKIYACEVRRKHHLERKMGGGFDSTEYDQQPQKYVSTFSKKIAPYASVIINGIYWAVNSPKLLTIPDAKYLLQPANTPWLPSSVGSPALPHRMLAICDISADPGGSIEFMNECTTIDTPFCLYDANRNKDTKSFSGPGVLVCSIDNMPTQLPLEATDFFGELLLPYLPDILQSDATKPFEEHNFSKTVHGSVIASNGKLTPNFKYITELR-NASKLRSKAQQMAASHDRKVLVLGAGYVSGPLVEYLTREGGNVSVTVASALKEEADGLASRFPGVEPVLLDVQENPDQLGKLIKSANVVVSLLPYYLHRSIAEQCVEHKTNMVTASYCVPALAELHQGALDADITCVNEVGLDPGIDHLLAMECFDDVHSSGGKVESFVSYCGGLPAPECSDNPLRYKFSWSPRGVLGNVLGSAKYLRNGKVVDIPAGG-LLDSRKSLDFLPGFALESLPNRDSTPYREKYGILEAHTVERHSLRYHGFCEALKGLVQIGMIDSNPHPSLHPMGPEITWRNLICNMLGLPDSNIFYENLKVEVLKRVEDNLKRMEAIEALGLLSDDKVVLCKTPLDSLSLYLSKRLALGEGERDVVIMRHEVDILWPDKRKEKRLINFVTYGDPKGYSAMAKTVGYPTAIATKMLLDGEIQTKGMVLPFTPEIYRPMLSRLRSEGISAVEKSIWI
#411976	Acerentomon_sp_9999213752	0.0	58.69	915	LDSVRGKTLAIRREEQSIWERRAPLSPTHVRKLVKAGVRVIVQPSNRRAYPMQTYANAGAIVKEDISEAPVVFGVKQVPVDALQPNKTYCFFSHTIKAQEGNMPLLDAILKKNIRLIDYEKMVVQSGQRIVAFGKYAGVAGMINILHGLGLRLLALGHHTPFMHIGPAHNYRNSMMARQAVRDAGYQIALGLLPRSIGPLTIVFTGTGNVSQGAQEIFQELPHEYIPVEMLQKVAEHGATNKIYACEVRRKHHLERKMGGGFDSTEYDQQPQKYVSTFSKKIAPYASVIINGIYWAVNSPKLLTIPDAKYLLQPANTPWLPSSVGSPALPHRMLAICDISADPGGSIEFMNECTTIDTPFCLYDANRNKDTKSFSGPGVLVCSIDNMPTQLPLEATDFFGELLLPYLPDILQSDATKPFEEHNFSKTVHGSVIASNGKLTPNFKYITELRNASKLRSKAQQMAASHDRKVLVLGAGYVSGPLVEYLTREGGNVSVTVASALKEEADGLASRFPGVEPVLLDVQENPDQLGKLIKSANVVVSLLPYYLHRSIAEQCVEHKTNMVTASYCVPALAELHQGALDADITCVNEVGLDPGIDHLLAMECFDDVHSSGGKVESFV---SYCGGLPAPECSDNPLRYKFSWSPRGVLGNVLGSAKYLRNGKVVDIPA-GGLLDSRKSLDFLPGFALESLPNRDSTPYREKYGILEAHTVERHSLRYHGFCEALKGLVQIGMIDSNPHPSLHPMGPEITWRNLICNMLGLPDSNIFYENLKVEVLKRVEDNLKRMEAIEALGLLSDDKVVLCKTPLDSLSLYLSKRLALGEGERDVVIMRHEVDILWPDKRKEKRLINFVTYGDPKGYSAMAKTVGYPTAIATKMLLDGEIQTKGMVLPFTPEIYRPMLSRLRSEGISAVEKS

# the subject IDs are made by script parse_ncbi_deflines_fasta, and are >[taxon]_[GI],
# although if no GI number parsed, meaningless transcript number is used.

# ATP6	Liophloeus_tessulatus_19877	4e-93	65.62	224	MTNLFSTFDPSTN-FNLSLNWLSIFLGIFIIPPMYWLIPSRMNMFWIKMIYTLHNEFKNLINKKEFKGSTLIFVSLFNLILFNNFLGLFPYIFTCTSHMTLTLSLSFPLWVTFMLNGWIKNTKFMFAHLVPQGAPSILLPFLVIIETISNIIRPGTLAIRLTANMITGHLLVTLLGNSGSSLNMYMLNILIFIQILLLILESAVAVIQSYVFSILATLYSSEVQ
# ATP6	Rhodopsona_rubiginosa_3200	4e-93	65.78	225	MMNNLFSIFDPSTNLMNIPFNWLSTFIGLMFIPFTFWLIPNRHFFFWNFILNKLHSEFKTLLGNINSNGSTFIFISLFSFILFNNFLGLFPYIFTSTSHLNLSLSISLPLWISFMLFGWINNSQHMFSHMIPQGTPNILMPFMVLIETISNIIRPSTLAVRLTANMIAGHLLMTLLSGTGNMIPSYLIMMLIIIQILLLILESAVAIIQSYVFAILSTLYSSEVN




# EOG502V7Z|DMELA_5.40|tsa.GAAB.1.fsa_nt.b.orf|Agrilus_planipennis_429235104|1|HNIQWIGAIRNISVLGRGSAFEVREESGPNLKSALRHILSVDLRDELIFPSCGTLFQLTTEFAGLGGNVGFLKNDAYLQTNYSIAEDFVIQTTFHGGYLRGLDNDLKVTLSDLFFLGGPLSLRGFHTRGIGPHAEGDALGATWFC

	if( $line =~ /^([^\|]+)\|([^\|]+)\|([^\|]+)\|([^\|]+)\|([^\|]+)\|([^\|]+)/ )
		{
		my $query_protein_ID = $1; my $subject = $4; my $seq = $6;my $length = length( $seq );


	#	print "\n$line\n";
	#	print "query_protein_ID:$query_protein_ID subject:$subject length:$length seq:$seq\n";		

		my $transcript_number;my $subject_GI;# subject gi is actually the species name


		#Occasjapyx_japonicus_999973103
		# more precisly, binomial followed by transcript number 9999\d+

		$blast_output_lines_parsed++;
		if($subject =~ /(.+)_([^_]+)$/)
			{
			$subject_GI = $1;$transcript_number = $2
			}else{
			die "\nerror 16, expecting /(.+)_([^_]+)/ \nsubject:($subject)\nline:$line\n"
			};
# print "subject:$subject subject_GI:$subject_GI transcript_number:$transcript_number\n";
		$record_all_subject_names{$subject_GI}=1;# subject_GI is actually taxonomjc name

		# if contains paralogs, 
		# for given query_ortholog-$subject_ome, will have multiple transcripts,
		# test that here.
		unless($test_paralogs{"$query_protein_ID.$subject_GI"} =~ /\t$transcript_number\t/)
			{$test_paralogs{"$query_protein_ID.$subject_GI"} .= "\t$transcript_number\t"};


		my $assign_code;
		if(exists($new_IDs{$subject_GI}))	# not used now
			{
			$assign_code = $new_IDs{$subject_GI};
			}else{
			$assign_code = "CODE" . $speciesnumber;$speciesnumber++;
			$new_IDs{$subject_GI} = $assign_code;
			};
	#	print "subject:$subject subject_GI:$subject_GI\n";		

		# phylotreepruner requires a uniqe ID for each member,
		# although currently i have one member per species,
		# so just create a unique code for each species.

		if($format == 2)
			{$subject_GI = $subject_GI . "|" . $assign_code}

		if($format == 3)
			{
			$subject_GI = $subject_GI . "|" . $transcript_number;
			}

		

		if($dont_print_duplicates == 1)
			{

			# $query_protein_ID is unique for each protein, different core species have same one. so ok to fileter.



			if(exists($subjects_printed{$query_protein_ID.$subject}))
			{
		#	print "already have PROT_ID:$query_protein_ID SUBJECT:$subject, subject_GI:$subject_GI\n";

			# i wasnt sure if blast orders hits best first, 
			# which would give sinmplere implementation.
			# i found cases it seems not to.
			if($subjects_printed{$query_protein_ID.$subject} > $length )
				{
			#	print "\tOLD is longer (prev:$subjects_printed{$query_protein_ID.$subject} new:$length)\n";
				}else{
			#	print "\tNEW is longer (prev:$subjects_printed{$query_protein_ID.$subject} new:$length)\n";
			#	$out_data{$query_protein_ID} .= ">$subject_GI\n$seq\n";

			#	print "1:out_data{query_protein_ID}:\n($out_data{$query_protein_ID})\n";
			# tried replacing old seq witrh new, but regexes i tried dont work, perhaps cause of | in name?
			#	# 2016-08-15: replace old seq, found longer one
			#	if($out_data{$query_protein_ID} =~ s/\>$subject_GI[\n\r]+\S+/">$subject_GI\n$seq\n"/se)
			#		{
			#		}else{die "\nerror trying to replace old short sequence for protein $query_protein_ID\n"};
			#	print "2:out_data{query_protein_ID}:\n($out_data{$query_protein_ID})\n";

				$subjects_printed{$query_protein_ID.$subject} = $length;
				$subjects_printed2{$query_protein_ID.$subject} = ">$subject_GI\n$seq\n";			

				};



			}else{
		#	print "dont yet have $query_protein_ID.$subject. storing for out data and recording length ($length)\n";
		#	$out_data{$query_protein_ID} .= ">$subject_GI\n$seq\n";
			$subjects_printed{"$query_protein_ID.$subject"} = $length;			
			$subjects_printed2{"$query_protein_ID.$subject"} = ">$subject_GI\n$seq\n";			
				# key looks like ATP6Rhodopsona_rubiginosa_3200


			};



			}else{#if($dont_print_duplicates == 1)
			die "\n\$dont_print_duplicates == 0 not currently supported. see notes in script...\n";
			$out_data{$query_protein_ID} .= ">$subject_GI\n$seq\n";
			}



		}else{# if line =~ blah blah

		if($line =~ /\w/){die "\nparse error:$line\n"};
		
		};
	#print $line;
	}
close IN;

print "
read blastP out file ($in)
	blast_output_lines_parsed:$blast_output_lines_parsed
";



if($dont_print_duplicates == 1)
	{
	my @best_seqs_keys = keys %subjects_printed2;# # key looks like ATP6Rhodopsona_rubiginosa_3200
	foreach my $index (0 .. $#best_seqs_keys)
		{
		my $key = $best_seqs_keys[$index];
	#	print "key:$key\n";

		# 2017 may, previiusyl had this:
	#	if($key =~ /^(\d+)\w/)
		# but not working in new proejct , 
		# seems changed from nunber to abreviated gene name
		if($key =~ /^([^\.]+)/)
			{
			my $query_protein_ID = $1;
			my $entry = $subjects_printed2{$key};
			$out_data{$query_protein_ID} .= $entry;
			}else{
			die "\nerror, ($index of $#best_seqs_keys), expecting query protein id to be \d+, then binomial, got:$key\n";
			};
		};
	};


##############

my @query_orth_subject_omes = keys %test_paralogs; @query_orth_subject_omes = sort @query_orth_subject_omes;
foreach my $qoso(@query_orth_subject_omes)
	{
	# $test_paralogs{"$query_protein_ID.$subject_GI"} =~ /\t$transcript_number\t/)
	my $hits = $test_paralogs{$qoso};#print "query_orth_subject_omes combination:$qoso hits:$hits\n";
	my $sp;my $orth;
	if($qoso =~ /(.+)\.(.+)/){$orth = $1; $sp = $2}else{die "\nerror 297\n"};

	if($hits =~ /\d\t\t\d/)
		{
		$species_for_which_paralogs_found{$sp} = 1;
		$orthologs_for_which_paralogs_found{$orth} = 1;
		# print "multiple transcripts found\n"; # die "";
		};

	};

my @paralog_species = %species_for_which_paralogs_found;
my @paralog_orths = %orthologs_for_which_paralogs_found;
print "number of species with paralogs:$#paralog_species, @paralog_species\n";
print "number of orthologs with paralogs:$#paralog_orths\n";


##############





my @keys = keys %out_data;#print "\n\@keys:@keys\n";die;
# keys array should look like:
#  100791 100111 100344 100568 

@keys = sort @keys;
my @all_subject_taxa = keys %record_all_subject_names;
@all_subject_taxa = sort @all_subject_taxa;

open(TABLE , ">parse_ortholog_results_TABLE") || die "";
print TABLE "Ortholog\tAnnotation\t";
foreach my $name(@all_subject_taxa)
	{
	print TABLE "$name\t";
	};
print TABLE "\n";

foreach my $ortholog( @keys )
	{
	my $ortholog_data = $out_data{$ortholog};
	my $ortholog_data_COPY  = $ortholog_data;
	my $ortholog_data_COPY2  = $ortholog_data;

	#print "\n\tortholog:$ortholog\northolog_data:$ortholog_data\n";	
	my %tax_for_this_otholog = ();
	while($ortholog_data_COPY =~ s/>([A-Z][a-z]+_[a-z]+)//)
		{my $taxon = $1;$tax_for_this_otholog{$taxon}++};
	my @unique_sp_for_this_ortholog = keys %tax_for_this_otholog;
	@unique_sp_for_this_ortholog = sort @unique_sp_for_this_ortholog;
	my $count_tax_this_orth = scalar @unique_sp_for_this_ortholog;
	#print "count_tax_this_orth:$count_tax_this_orth SP:@unique_sp_for_this_ortholog\n";

	my %IDs_for_this_otholog = ();
	while($ortholog_data_COPY2 =~ s/>([A-Z][a-z]+_[a-z]+.+)//)
		{
		my $ID = $1;$IDs_for_this_otholog{$ID}++;#print "ID:$ID\n";
		if($IDs_for_this_otholog{$ID} >= 2){die "\nerror 277 ... ortholog:$ortholog has repeated transcript $ID\n"};
		};



	if($count_tax_this_orth >= $ortholog_taxon_limit )
		{

	if($ortholog =~ /^[WY]P_\d+/)
		{
		open(OUT, ">$string.$ortholog") || die "error cant open out for ortholog ($ortholog)\n";
		print OUT $ortholog_data;
		close OUT;
		}elsif($ortholog =~ s/^(\d+)$/WP_$1/)
		{
		open(OUT, ">$string.$ortholog") || die "error cant open out for ortholog ($ortholog)\n";
		print OUT $ortholog_data;
		close OUT;
		}elsif($ortholog =~ /^EOG[A-Z0-9]+$/)
		{
		# EOG5V41PN, 1KITE's
 
		open(OUT, ">$string.$ortholog") || die "error cant open out for ortholog ($ortholog)\n";
		print OUT $ortholog_data;
		close OUT;
		

		}else{
		print  "42. $ortholog != /[WY]P_\d+\n";
		open(OUT, ">$string.$ortholog") || die "error cant open out for ortholog ($ortholog)\n";
		print OUT $ortholog_data;
		close OUT;
		
		}

		}else{
		$filter_orthologs_due_to_not_many_tax++;
		}

	$ortholog =~ s/^WP_//;
	my $annotation = "NA";
	if($store_annotations{$ortholog} =~ /\w/)
		{
		$annotation = $store_annotations{$ortholog};
		};

	print TABLE "$ortholog\t$annotation\t";


	foreach my $name(@all_subject_taxa)
		{
		if(exists($tax_for_this_otholog{$name}))
			{print TABLE "1\t" 
			}else{print TABLE "0\t"};
		}
	print TABLE "\n";

	};#foreach my $ortholog(@keys)




print "
count orthologs:$#keys
each printed to its own file
filter_orthologs_due_to_not_many_tax:$filter_orthologs_due_to_not_many_tax
 (ortholog_taxon_limit:$ortholog_taxon_limit)



FIN.
";








############################################################################################################

sub read_command_arguments
{

my $arguments = shift;
$arguments =~ s/\n//;$arguments =~ s/\r//;


# $dont_print_duplicates =1; $format = 3
# perl ~/usr_scripts/parse_ortholog_results.pl -in insecta_ortholog_hits.all -out_prefix insectaNUCL

# $in = $ARGV[0];
# $string = $ARGV[1];


if($arguments =~ /-in\s+(\S+)/)
	{
	$in = $1;
	open(INTEST , $in) || die "\ninput not found in current directlry:$in\n";
	close INTEST;
	}else{
	print "\nerror reading command arguments  (-in)\n";die""
	}

if($arguments =~ /-out_prefix\s+(\S+)/)
	{
	$string = $1;
	}else{
	print "\nNO out prefix?\n";die;
	};

if($arguments =~ /-filter_duplicates/)
	{
	$dont_print_duplicates = 1;
	}


if($arguments =~ /-output_id_format\s+(\S+)/)
	{
	$format = $1;
	}else{
	die "";
	};


};

############################################################################################################



