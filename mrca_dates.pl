
$in = $ARGV[0]; 	# reference phylogeny from which to get dates
$in2 = $ARGV[1]; 	# list_constrained_IDs
$out = $ARGV[2];	# outfile name


# needed because oddly, Kawahara et al 2019 appear to divide time calibrated dates in MR versions by 100
$scale_dates = 0; # default 0
$scale_factor = 1;


###########################################################################################
open(IN2, $in2) || die "\nerror cant open file $in2\n";
while(my $line = <IN2>)
	{
	$line =~ s/\n//;$line =~ s/\r//;
	if($line =~ /\w/)
		{
		my $synthesis_sp = $line;
		if($synthesis_sp =~ /_sp$/)
			{}else{
			if($synthesis_sp =~ /^(\w+)_(\w+)$/)
				{
				my $genus = $1; my $species = $2;$exemplar_sp_of_genus{$genus} = $species;$synthesis_exemplars_stored++;
				$synthesis_species{$synthesis_sp} = 1;
				};
			};
		$all_synthesis_tree_species{$line} = 1; # print "line:$line\n";
		};
	
	};

print "
synthesis_exemplars_stored:$synthesis_exemplars_stored
";
###########################################################################################
# read time calibrated backbone from file
my $tree;
$trees_parsed =0;
open(IN, $in) || die "\nerror cant open file $in\n";
while(my $line = <IN>)
	{
	$line =~ s/\n//;$line =~ s/\r//;
	if($line =~ /\(/){$tree = $line;$trees_parsed++};
	};
close IN;
unless($trees_parsed == 1){die "\nerror unexpected number of newick strings\n"}
$treelength = length($tree);
print "\ntree of length $treelength read from file $in\n";
###########################################################################################
# parse newick

# print "tree:$tree\n";die "";
unless($tree =~ /\:\d/){die "\nreference phylogeny should contain branchlengths. quitting.\n"};
$node =0;
while($tree =~ /[\(]([\w\_]+)\:([\d\.]+)\,([\w\_]+)\:([\d\.]+)\)\:([\d\.]+)/)
	{
	my $child1_nodeID = $1; my $child1_branchlength = $2;	my $child2_nodeID = $3; my $child2_branchlength = $4; my $adjoining_branchlength = $5;
#	print "child1_nodeID:$child1_nodeID child1_branchlength:$child1_branchlength child2_nodeID:$child2_nodeID child2_branchlength:$child2_branchlength adjoining_branchlength:$adjoining_branchlength\n";

	# branchlength of new node will be addition of 2:
	my $new_length = $child1_branchlength+$adjoining_branchlength;

	my $ref_terminal_ID_1 = $child1_nodeID;
	if($child1_nodeID =~ /^(\w+)_(\w+)$/)
		{
		my $genus = $1; my $species = $2;
		if($species eq "sp" && $exemplar_sp_of_genus{$genus} =~ /\w/){$ref_terminal_ID_1 = $exemplar_sp_of_genus{$genus};$exemplar_switches++};
		};
	my $ref_terminal_ID_2 = $child2_nodeID;
	if($child2_nodeID =~ /^(\w+)_(\w+)$/)
		{
		my $genus = $1; my $species = $2;
		if($species eq "sp" && $exemplar_sp_of_genus{$genus} =~ /\w/){$ref_terminal_ID_2 = $exemplar_sp_of_genus{$genus};$exemplar_switches++};
		};

	if( $mrcas{$child1_nodeID . "\t" . $child2_nodeID} =~ /\d/)
		{
		die "\nerror 34.\n"
		}else{
		};
	$mrcas{$ref_terminal_ID_1 . "\t" . $ref_terminal_ID_2} = $child1_branchlength;

	
	# note: replacing internal node with a randomly selected terminal, not a new internal node id as often done elsewhere.
	# this is since we need to refer to terminals derived from internal nodes.

	my $prefered_terminal = $child2_nodeID;
	if($synthesis_species{$child1_nodeID} == 1)
		{
		$prefered_terminal = $child1_nodeID;
		}
#	$tree =~ s/[\(]([\w\_]+)\:([\d\.]+)\,([\w\_]+)\:([\d\.]+)\)\:([\d\.]+)/$child1_nodeID:$new_length/;
#	$tree =~ s/[\(]([\w\_]+)\:([\d\.]+)\,([\w\_]+)\:([\d\.]+)\)\:([\d\.]+)/$child2_nodeID:$new_length/;
	$tree =~ s/[\(]([\w\_]+)\:([\d\.]+)\,([\w\_]+)\:([\d\.]+)\)\:([\d\.]+)/$prefered_terminal:$new_length/;

	$node++;
	};

print "\nremaining newick string:$tree\n";
print "
exemplar_switches:$exemplar_switches
";
my $parsed_newick_lenght = length($tree);
if($parsed_newick_lenght >= 100)
	{die "\nerror, newick string not parsed, this can happen if your time calibrated backbone has multifurcations, please removed these.\n"}
###########################################################################################

# list mrcas

open(OUT, ">$out") || die "";
open(OUT2, ">$out.reduced") || die "";
$calibration_number = 0;
my @all_mrcas = keys %mrcas; @all_mrcas = sort @all_mrcas;
foreach my $node(@all_mrcas)
	{
	my @both = split /\t/ , $node;
	my $date = $mrcas{$node};
#	print "mrca of $node is dated $date\n";

	if($scale_dates == 1){ $date = $date * $scale_factor};


	if($all_synthesis_tree_species{$both[0]} == 1 && $all_synthesis_tree_species{$both[1]} == 1)
		{
	#	print "\tCAN be applied to synth tree\n";
		$applicable++;
		#mrca = ANTHOPHILA Philanthus_gibbosus Bombus_impatiens
		#min = ANTHOPHILA 1
		#max = ANTHOPHILA 1
		print OUT "mrca = CALIB$calibration_number $both[0] $both[1]\n";
		print OUT "min = CALIB$calibration_number $date\n";
		print OUT "max = CALIB$calibration_number $date\n";
		$calibration_number++;

		if($donealready{$both[0]} == 1 || $donealready{$both[1]} == 1)
			{}else{
			print OUT2 "mrca = CALIB$calibration_number $both[0] $both[1]\n";
			print OUT2 "min = CALIB$calibration_number $date\n";
			print OUT2 "max = CALIB$calibration_number $date\n";
			};
		$donealready{$both[0]} = 1;$donealready{$both[1]} = 1;

		}else{
	#	print "\tCAN NOT be applied to synth tree\n";
		$not_applicable++;
		
		};
	};

print "
applicable:$applicable
not_applicable:$not_applicable
";

if($applicable <= 5)
	{
	print "\nnot many dates parsed, you might try alternative option on line 59\n";
	};

close OUT;
close OUT2;

###########################################################################################















