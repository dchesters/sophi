
# $database_file 	= $ARGV[0];
# $key_file 	= $ARGV[1];


###################################################################################################################################
#
#
#
# 	parse_taxon_from_fastafile.pl, 
#	Perl script for consolidating two specific files, a fasta database and a file of taxonomic information, 
#		taxa in the latter are parsed from the db and printed to outfile
#
#    	Copyright (C) 2013-2016  Douglas Chesters
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
# 	CHANGE LOG
# 	jan2010: print LOG
# 	june 2011: slighly different format used for input sequence fastafile.
# 	april 2013: dont read whole fasta file into memory. no good for larger databases.
# 	apr2013: make a new key file removing species not included in database.
# 	oct2013: input files given in command line. include GPL.
# 	10oct2013: option for family name on output IDs
#	25aug2014: option to ignore unlabelled sp
#	31dec2014: you can limit to a list of species, given in array below
#		used to get genomic sequecnes only
#	
#	
#	
#	
#	
#	
#
#
# 	NOTES
# 	blastclust -i clade.IA.fas -o clade.IA.blastclust_outfile.30 -S 30 -p F -e F -a 2
# 	wc -l inv.fas #4,450,986
# 	split -l 1000000 inv.fas
# 	4,450,986
# 	cat xab.parsed >> xaa.parsed
#
#
#
#
#
#	
#	
#	



my $arg_string  = join ' ', @ARGV;
#####################################
read_command_arguments($arg_string);#
#####################################



$parsed_database_name 	= "$database_file.parsed";

# $memory_efficient	= 1; 	# 1= read only one entry into memory at a time, allows large fasta databases to be used.

# $id_format 		= 3; 	# 1 = tobycode, underscore, accession. 
				# 2 = ncbi_taxon_number, underscore, accession.
				# 3= full species name, underscore, accession
				# 4= family, underscore, full species name, underscore, accession

# $parse_binomial_labelled_only =1;# default = 0 

	# if this is empty, proceeds as normal. 
	# if you write some species names in the array, seqeunces only for these will be parsed




@parse_these_species_only = (
#"Acanthocasuarina_muellerianae","Acanthosoma_haemorrhoidale","Acromyrmex_echinatior","Acyrthosiphon_pisum","Adineta_vaga","Agrilus_planipennis","Aleochara_curtula","Amphimedon_queenslandica","Anopheles_albimanus",
#"Anopheles_arabiensis","Anopheles_atroparvus","Anopheles_christyi","Anopheles_culicifacies","Anopheles_dirus","Anopheles_epiroticus","Anopheles_farauti","Anopheles_funestus","Anopheles_maculatus","Anopheles_melas","Anopheles_merus",
#"Anopheles_minimus","Anopheles_quadriannulatus","Anopheles_sinensis","Anopheles_stephensi","Anoplophora_glabripennis","Anurida_maritima","Apachyus_charteceus","Apis_dorsata","Apis_florea","Apis_mellifera","Aplysia_californica",
#"Aposthonia_japonica","Aretaon_asperrimus","Ascaris_suum","Atelura_formicaria","Athalia_rosae","Atta_cephalotes","Bactrocera_cucurbitae","Bibio_marci","Biomphalaria_glabrata","Bittacus_pilicornis","Blaberus_atropos",
#"Blattella_germanica","Bombus_impatiens","Bombus_terrestris","Bombylius_major","Boreus_hyemalis","Botryllus_schlosseri","Branchiostoma_floridae","Brugia_malayi","Caenorhabditis_angaria","Caenorhabditis_brenneri","Caenorhabditis_briggsae",
#"Caenorhabditis_elegans","Caenorhabditis_japonica","Calopteryx_splendens","Campodea_augens","Camponotus_floridanus","Capitella_teleta","Capsaspora_owczarzaki","Centruroides_exilicauda","Cephus_cinctus","Cerapachys_biroi",
#"Ceratitis_capitata","Ceratophyllus_gallinae","Ceratosolen_solmsi","Cercopis_vulnerata","Chironomus_tentans","Chrysis_viridula","Cimex_lectularius","Ciona_intestinalis","Clonorchis_sinensis","Conwentzia_psociformis",
#"Copidosoma_floridanum","Cordulegaster_boltonii","Corydalus_cornutus","Cosmioperla_kuna","Cotesia_vestalis","Crassostrea_gigas","Cryptocercus_wrighti","Ctenocephalides_felis","Danaus_plexippus","Daphnia_pulex",
#"Dendroctonus_ponderosae","Dermatophagoides_farinae","Diaphorina_citri","Drosophila_albomicans","Drosophila_biarmipes","Drosophila_bipectinata","Drosophila_elegans","Drosophila_eugracilis","Drosophila_ficusphila","Drosophila_kikkawai",
#"Drosophila_melanogaster","Drosophila_miranda","Drosophila_pseudoobscura","Drosophila_rhopaloa","Drosophila_simulans","Drosophila_suzukii","Drosophila_takahashii","Drosophila_willistoni","Drosophila_yakuba","Dyseriocrania_subpurpurella",
#"Echinococcus_granulosus","Echinococcus_multilocularis","Ectopsocus_briggsi","Elaeophora_elaphi","Empusa_pennata","Ephemera_danica","Ephemera_danica","Epiophlebia_superstes","Essigella_californica","Euroleon_nostras",
#"Eurytemora_affinis","Folsomia_candida","Fonticula_alba","Fopius_arisanus","Forficula_auricularia","Frankliniella_occidentalis","Frankliniella_cephalica","Galloisiana_yuasai","Globodera_pallida","Glossina_austeni",
#"Glossina_brevipalpis","Glossina_fuscipes","Glossina_pallidipes","Grylloblatta_bifratrilecta","Gynaikothrips_ficorum","Gyrinus_marinus","Gyrodactylus_salaris","Haemonchus_contortus","Halyomorpha_halys","Haploembia_palaui",
#"Harpegnathos_saltator","Heliconius_melpomene","Helobdella_robusta","Heterorhabditis_bacteriophora","Homalodisca_vitripennis","Hyalella_azteca","Hydra_magnipapillata","Hymenolepis_microstoma","Inocellia_crassicornis",
#"Isonychia_bicolor","Ladona_fulva","Lasioglossum_albipes","Latrodectus_hesperus","Leptinotarsa_decemlineata","Leptopilina_clavipes","Limnephilus_lunatus","Limulus_polyphemus","Linepithema_humile","Lipara_lucens",
#"Liposcelis_bostrychophila","Loa_loa","Lottia_gigantea","Lucilia_cuprina","Lutzomyia_longipalpis","Lytechinus_variegatus","Machilis_hrabei","Manduca_sexta","Mantis_religiosa","Mastotermes_darwiniensis",
#"Mayetiola_destructor","Megachile_rotundata","Meinertellus_cundinamarcensis","Melitaea_cinxia","Meloe_violaceus","Menopon_gallinae","Metallyticus_splendidus","Metaseiulus_occidentalis","Microplitis_demolitor","Mnemiopsis_leidyi",
#"Musca_domestica","Nasonia_giraulti","Nasonia_longicornis","Nasonia_vitripennis","Necator_americanus","Nemophora_degeerella","Nilaparvata_lugens","Nilaparvata_lugens","Notostira_elongata","Occasjapyx_japonicus",
#"Okanagana_villosa","Onchocerca_volvulus","Oncopeltus_fasciatus","Onthophagus_taurus","Opisthorchis_viverrini","Orussus_abietinus","Orussus_abietinus","Osmylus_fulvicephalus","Pachypsylla_venusta","Panagrellus_redivivus",
#"Panorpa_vulgaris","Parasteatoda_tepidariorum","Parides_eurimedes","Patiria_miniata","Periplaneta_americana","Perla_marginata","Peruphasma_schultei","Phlebotomus_papatasi","Planococcus_citri","Platycentropus_radiatus",
#"Plutella_xylostella","Pogonomyrmex_barbatus","Polyommatus_icarus","Priapulus_caudatus","Prorhinotermes_simplex","Prosarthria_teretrirostris","Pseudomallada_prasinus","Ranatra_linearis","Rhodnius_prolixus","Rhyacophila_fasciata",
#"Saccoglossus_kowalevskii","Schistosoma_haematobium","Schistosoma_mansoni","Schmidtea_mediterranea","Sminthurus_viridis1","Solenopsis_invicta","Stegodyphus_mimosarum","Steinernema_carpocapsae","Steinernema_feltiae",
#"Steinernema_glaseri","Steinernema_monticolum","Steinernema_scapterisci","Stenobothrus_lineatus","Strigamia_maritima","Strongylocentrotus_purpuratus","Stylops_melittae","Tenthredo_koehleri","Tetranychus_urticae","Tetrix_subulata",
#"Tetrodontophora_bielanensis","Thermobia_domestica","Thrips_palmi","Timema_cristinae","Trialeurodes_vaporariorum","Triarthria_setipennis","Trichinella_spiralis","Trichocera_saltator","Trichogramma_pretiosum","Tricholepidion_gertschi",
#"Trichuris_muris","Trichuris_suis","Trichuris_trichiura","Triodia_sylvina","Velia_caprai","Xanthostigma_xanthostigma","Xenophysella_greensladeae","Yponomeuta_evonymella","Zootermopsis_nevadensis","Zorotypus_caudelli","Zygaena_fausta"
);



$logfile 	= "parse_order_from_endop_fastafile_LOG";







open(LOG , ">$logfile") || die "error 31. cant open file\n";
open (IN, "$key_file") || die "error 32. cant open file:$key_file\n";

%storelabels = ();
$taxid;

my $count_keys = 0;
while (my $line= <IN>)
	{
	# LZ1M4Mi7cal Micropterix_calthella 41027 species suborder:Zeugloptera family:Micropterigidae genus:Micropterix species:Micropterix calthella 

	if($line =~ /^(\S+)\s(\S+)\s(\d+)\s(.+)/)
		{
		my $tobycode = $1;my $full_species_name = $2;my $taxid = $3;my $rest = $4;
	#	if($line =~ /1006788/){print "\n\n(taxid:$taxid) line in keyfile $line\n"}

		my $proceed = 1;
		if($parse_binomial_labelled_only == 1)
			{
			if($full_species_name =~ /^[A-Z][a-z]+_[a-z]+$/)
				{
				if($full_species_name =~ /^[A-Z][a-z]+_sp$/ || $full_species_name =~ /^[A-Z][a-z]+_spp$/ || $full_species_name =~ /^[A-Z][a-z]+_cf$/ || $full_species_name =~ /^[A-Z][a-z]+_nr$/ )
					{$proceed = 0}else{}
				}else{$proceed = 0}
			};

		if(scalar @parse_these_species_only >=2)
			{
			$proceed = 0;
			foreach my $species98(@parse_these_species_only)
				{
				if($full_species_name eq $species98){$proceed = 1;print "genomic:$species98\n"}
				};
			};

		if($proceed == 1)
			{

		my $family_name = "NA";if($rest =~ / family:(\w+) /){$family_name = $1}

		if($id_format =~ /[34]/)
			{
			$storelabels{$taxid} = $full_species_name;			
			}else{
			$storelabels{$taxid} = $tobycode;# ncbi tax number is key, tobycode or species name is assigned.
			}
		$storefam{$taxid} = $family_name;
		$count_keys++;

		if($line =~ /^(\S+)\s\S+\s\d+\sspecies\s/)
			{$species_level_codes{$tobycode} = 1}

			};

		}else{
		#if($line =~ /^(\S+)/){die "BUG1\n$line\n"}
		}
	};

close(IN);


my @species_keys = keys %species_level_codes;


print "\nkey file $key_file has been read and contains $count_keys taxa, $#species_keys species level tobycodes\n";
print LOG "\nkey file $key_file has been read and contains $count_keys taxa, $#species_keys species level tobycodes\n";



unless ($#species_keys >= 1){die "\nerror 194\n"};
unless ($count_keys >= 1){die "\nerror 195\n"};



my $count_matching_taxids = 0;

my $total_sequences_in_DB =0;

	open(LOGMISSING, ">missingtaxids") || die "\nerror 78\n";

%fasta_seqs;
$length_NOT_included =0;$length_included = 0;




if($memory_efficient == 1)
	{
	open(OUT45, ">$parsed_database_name") || die "\nerror\n";
	close(OUT45);
	}

#########################
fasta_format_to_nexus();#
#########################

if($memory_efficient == 1)
	{}


if($memory_efficient == 0)
	{
###################
print_parsed_db();#
###################
	}

#@species_keys = keys %species_level_codes;


foreach my $sp(@species_keys)
	{
	if($species_level_codes{$sp} == 1)
		{
		$count_species_with_sequence_NOT_present_in_db++;
		print LOG "species $sp not present in DB\n";
		}
	if($species_level_codes{$sp} == 2)
		{
		$count_species_with_sequence_present_in_db++;
	
		}
	}

print "count_species_with_sequence_NOT_present_in_db ($count_species_with_sequence_NOT_present_in_db)
count_sequences_printed_to_output:$count_sequences_printed_to_output
";

unless($count_sequences_printed_to_output >= 1){die "\n\nERROR no seqs printed\n\n"};

print LOG "count_species_with_sequence_NOT_present_in_db ($count_species_with_sequence_NOT_present_in_db)
count_species_with_sequence_present_in_db ($count_species_with_sequence_present_in_db)
";

$length_NOT_included =~ s/(\d\d\d)(\d\d\d)$/,$1,$2/;$length_included  =~ s/(\d\d\d)(\d\d\d)$/,$1,$2/;

print "length of data included in output ($length_included) bytes.\nlength of data ommited from output ($length_NOT_included) bytes. ";
print "these are sequences in the file $database_file which have ncbi taxon numbers not listed in $key_file. ";
print " ommited sequences are listed in file missingtaxids\n";
print LOG "length of data included in output ($length_included) bytes.\nlength of data ommited from output ($length_NOT_included) bytes. ";
print LOG "these are sequences in the file $database_file which have ncbi taxon numbers not listed in $key_file. ";
print LOG " ommited sequences are listed in file missingtaxids\n";

 
write_reduced_keyfile();



exit;



sub fasta_format_to_nexus
	{
	print "sub fasta_format_to_nexus. reading $database_file. this may take a minute\n";
	print LOG "sub fasta_format_to_nexus. reading $database_file.\n";

	my $count_missing_taxids =0;
	my $current_seq_length;
	my $missing_data_character = "N";
	my $count_no_fasta_entries = 0;

	open(FASTA_IN, $database_file) || die "Cant open input:$database_file.\n";
	my $file_as_string = ""; my @all_lines = ();

	my $fasta_entry = "";
	while(my $fasta_line= <FASTA_IN>)
		{
		#$file_as_string .= $line

		if($fasta_line =~ /^>(.+)/)
			{
			my $fasta_id = $1;$count_no_fasta_entries++;
			#print "$fasta_line";

			unless(length($fasta_entry)<=2)
				{
				$fasta_entry =~ s/^>//;
				########################################
				process_this_fasta_entry($fasta_entry);#
				########################################
				
				}
			$fasta_entry = $fasta_line;


			my $IDcopy = "NA"; if($fasta_id =~ /^([\w\d\_\-]+)/){$IDcopy = $1};
			if($count_no_fasta_entries =~ /00000$/)
				{print "entry ($IDcopy), total read in db:$count_no_fasta_entries, id's found:$count_matching_taxids\n";
				print "\tcount_matching_taxids:$count_matching_taxids\n"}

			}else{
			$fasta_entry .= $fasta_line;
			
			}




		};
	close(FASTA_IN);





####################################################################################################################

sub process_this_fasta_entry
{
my $line = shift;

if(length($line)<=30)
	{
	print LOG "why so short:$line\n .... ignoring ... \n"
	}else{
	#	if($line =~ /1006788/){print "A line in fasta $line\n"}

	if ($line =~ /^([^_]+)_(.+)\n/ )
		{
		my $current_accession = $1;my $current_taxid=$2;#$locus = $3;
		#print "current_accession:$current_accession current_taxid:$current_taxid\n";
		#current_accession:Strongylus current_taxid:vulgaris_689851437
		if($current_accession =~ /^[A-Z][a-z]+$/ && $current_taxid =~ /^[a-z]+_[0-9]+$/)		
			{
			die "\nerror 358
seems the file you are running this on, has IDs printed incorrectly. expecting >[accession]_[tax_ID]
probably need alteration to script create_fasta_database_from_genbank_flatfiles.pl, \$output_ID_format = 1
\n"
			};
		#	if($line =~ /1006788/){print "(cti:$current_taxid) B line in fasta $line\n"}

		#	print "checking if id $current_taxid is stored\n";
		my $entry_length = length($line);
		$line =~ s/^.+\n//;	$line =~ s/\n//;$line =~ s/\r//;$current_seq_length = length($line);
		$total_sequences_in_DB++;

		if(length($storelabels{$current_taxid})>=1)
			{
		#	if($line =~ /1006788/){print "C line in fasta $line\n"}

			$key = $storelabels{$current_taxid} . "_" . $current_accession;# . "__$locus";
			$count_matching_taxids++;

			if($memory_efficient == 0)
				{
				$fasta_seqs{$key} = $line;
				}else{
				open(OUT45, ">>$parsed_database_name") || die "\nerror\n";
				if($id_format =~ /[13]/)
					{print OUT45 ">$key\n$line\n"};
				if($id_format == 4){print OUT45 ">" , $storefam{$current_taxid} , "_$key\n$line\n"}
				if($id_format == 2)
					{print OUT45 ">$taxid" , "_$current_accession\n$line\n"}
				close(OUT45);
				$count_sequences_printed_to_output++;
				}
			$species_level_codes{$storelabels{$taxid}}  = 2;

			$length_included += $entry_length;
			}else{
			$count_missing_taxids++;print LOGMISSING "$current_accession $taxid\n";
			# do nothing, you proably dont want this taxon
			#print "$taxid doesnt exist\n";
			$length_NOT_included += $entry_length;

			}


			}else{
			print "ignoring entry:$line\n";
			#die "BUG2\nhash out this line if not many\n"
			};

		}


}



####################################################################################################################



#		}# for my $each_line(1 .. $#all_lines)


	};# sub fasta_format_to_nexus





close(LOGMISSING);

sub print_parsed_db
{
print "\nsub print_parsed_db\n";
my $noininput = scalar @all_lines;
	my @fasta_seqs_keys = keys %fasta_seqs;
	my $number_of_taxa = scalar @fasta_seqs_keys;
	print  "$number_of_taxa entries read into memory.\nin input:$noininput\n";
	print  LOG "$number_of_taxa entries read into memory.\n";
	print "count_missing_taxids:$count_missing_taxids\n";

print "\nof the $total_sequences_in_DB total seqeunces in the database ($database_file),"; 
print	" $count_matching_taxids sequence entries matched one of taxa in your keyfile ($key_file), ";
print "which presumably are a subset of the whole release division\n\n";

print LOG "\nof the $total_sequences_in_DB total seqeunces in the database ($database_file),"; 
print LOG " $count_matching_taxids sequence entries matched one of taxa in your keyfile ($key_file), ";
print LOG "which presumably are a subset of the whole release division\n\n";


	############################################
	@fasta_seqs_keys = sort @fasta_seqs_keys;
###########################################


	open(NEXUS_OUT ,">$parsed_database_name") || die "cant open output file:$format_conversion_output_file\n";
	
	for $i(0 .. $#fasta_seqs_keys)
		{
		$current_name = $fasta_seqs_keys[$i];
		$current_seq = $fasta_seqs{$current_name};
	#	$current_name =~ 
		# print "current_name:$current_name current_seq:$current_seq\n";
		if(length($current_name)<=1 || length($current_seq)<=1){print "warning: zero length of current entry. quitting\n";die}

		print NEXUS_OUT ">$current_name\n$current_seq\n";
		};

	close(NEXUS_OUT);
}





close(LOG);





###############################################################

sub write_reduced_keyfile
{

print "writing new key file ($key_file.reduced_according_to.$database_file) in which species that are not found in the database ($database_file) have been removed\n";


open (IN, "$key_file") || die "error 335. cant open file:$key_file\n";
open (OUT3, ">$key_file.reduced_according_to.$database_file") || die "error 336. cant open file\n";

my $count_keys = 0;
my $total_deleted =0;my $bold_deleted =0;
while ($line= <IN>)
	{
	# LZ1M4Mi7cal Micropterix_calthella 41027 species suborder:Zeugloptera family:Micropterigidae genus:Micropterix species:Micropterix calthella 

	if($line =~ /^(\S+)\s(\S+)\s(\d+)\s/)
		{
		$tobycode = $1;
		$taxid = $3;

		if($species_level_codes{$tobycode} == 1)
			{
		#	print "deleting $tobycode\n";
			$total_deleted++;if($tobycode =~ /BOLD/){$bold_deleted++}
			}else{
			print OUT3 $line
			}

		}
	}


close(IN);
close(OUT3);

my $nonBOLDrm = $total_deleted-$bold_deleted;
#print "total deleted ($total_deleted) bold deleted ($bold_deleted) nonBOLD removed ($nonBOLDrm) \n";

}


##############################################################################################################




sub read_command_arguments
{

my $arguments = shift;
$arguments =~ s/\n//;$arguments =~ s/\r//;

print "
";


if($arguments =~ /-fasta_db\s+(\S+)/)
	{
	$database_file = $1;
	}else{
	die "\n\nerror, please give name of fasta database as so: -fasta_db [file_name]\n\n";
	};
# $database_file 	= $ARGV[0];

# $key_file 	= $ARGV[1];
if($arguments =~ /-taxa_key\s+(\S+)/)
	{
	$key_file = $1;
	}else{
	die "\n\nerror, please give name of taxa key file as so: -taxa_key [file_name]\n\n";
	};

if($arguments =~ /-memory_efficient/)
	{
	$memory_efficient	= 1;
	}else{
	$memory_efficient	= 0;
	}
# $memory_efficient	= 1; 	# 1= read only one entry into memory at a time, allows large fasta databases to be used.

if($arguments =~ /-id_format\s+(\d)/)
	{
	$id_format = $1;
	}else{
	die "\n\nerror please give format for fasta IDs as so: -id_format [1/2/3/4]";
	};
# $id_format 		= 3; 	# 1 = tobycode, underscore, accession. 
				# 2 = ncbi_taxon_number, underscore, accession.
				# 3= full species name, underscore, accession
				# 4= family, underscore, full species name, underscore, accession

if($arguments =~ /-binomials_only/)
	{
 	$parse_binomial_labelled_only =1;# default = 0 
	}else{
 	$parse_binomial_labelled_only =0;# default = 0 
	};
# $parse_binomial_labelled_only =1;# default = 0 


# -fasta_db -taxa_key -id_format -memory_efficient -binomials_only


};


##############################################################################################################









