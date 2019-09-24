


# 
# 
# 
# change log 
# 2017 OCT 30: checks if pb outfile is parsed but no significance results found (though seems no cases of this)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 



$user_param = $ARGV[0];
$user_val = $ARGV[1];

unless($user_val =~ /\d/){die "\ncommadn errror\n"};

my @file_list = split /\n/ , `ls phylobayesOUT.*.clo_sample.ppht`;

my $count = scalar @file_list;
#print "@file_list";die;
unless($count >= 1){die "\nerror, cant find phylobayes results files (ls phylobayesOUT.*.clo_sample.ppht)\n"};

open(LOG , ">parse_phylobayes_LOG") || die "";

foreach my $file(@file_list)
	{
	my $pv = "NA"; my $zs = "NA";
	open(IN, $file) || die "";
	my $entry = "";
	while (my $line = <IN>)
		{
		$entry .= $line;
		#print "$line";
		if($line =~ /^p-value\s+\:\s+(.+)/)
			{
			$pv = $1;
			#p-value    : 1
			#z-score    : -2.81257
			};
		if($line =~ /^z-score\s+\:\s+(.+)/)
			{
			$zs = $1;
			};


# * Machilis_hrabei|709808518            0.02	2.329
		if($line =~ /\s\*\s+([A-Z][a-z]+_[a-z]+)/)
 			{my $tax  = $1;$whichtax{$tax}++}


		};
	close IN;


# high z-score = heterogeneity, low z-score = homogeneity, http://www.nature.com/articles/srep20528
# comparing list of values:high significance = high z-score, not signofocant = low (minus) zscore.


	my $passes=0;
	if($user_param eq "pv")
		{
		if($pv eq "NA")
			{
			$passes = 1;$no_test_score_found++;
			}
		elsif($pv <= $user_val)# test of compositional homogeneity
			{
			}else{
			$passes=1
			};		
		}elsif($user_param eq "zs")
		{
		if($zs eq "NA")
			{
			$passes = 1;$no_test_score_found++;
			}
		elsif($zs >= $user_val)
			{}else{$passes=1};		
		}else{
		die "\nerror 60.\n"
		};

my $filename = $file;
	if($passes == 0)
		{
		print "File:$file PV:$pv ZS:$zs\n";#$entry\n";
		$count_heterogeneous++;
		}else{
		$count_homogeneous++;
		print "File:$file PV:$pv ZS:$zs\n";#$entry\n";
	#File:phylobayesOUT.insectaNUCL.WP_646626466.clo_sample.ppht PV:1 ZS:-2.81257
		$string = "_pruned.fa"; 
		$file =~ s/phylobayesOUT\.(.+clo)_sample\.ppht/$1$string/;
		system("cp $file $file.RM1");
		}
	print LOG "File:$filename PV:$pv ZS:$zs passed_filter:$passes\n";
	};

my $perent = ($count_heterogeneous / ($count_heterogeneous+$count_homogeneous))*100;

$perent =~ s/(\d+\.\d)\d+/$1/;

print "\nnumber phylobayes results files:$count\n";

print "
	count_heterogeneous (removed):$count_heterogeneous
	count_homogeneous:$count_homogeneous
	no_test_score_found:$no_test_score_found
	percent removed:$perent
\n";
print LOG "\ncount_heterogeneous (removed):$count_heterogeneous
count_homogeneous:$count_homogeneous\n";

close LOG;

my @keys = keys %whichtax;
foreach my $key (@keys)
	{
	my $val = $whichtax{$key};
	#print "key:$key val:$val\n";
	};







