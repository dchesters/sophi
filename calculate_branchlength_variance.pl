



# 
# 
# 
# 
# writes temp files varianceR table_for_R R_Variance_OUT
# 
# 
# 
# 
# CHANGE LOG
# 
# 08 APR 2016: write R script anew in working folder, each time script invoked
# 23 JUN 2016: user specified output
# 17 AUG 2016: option to calculate av boot support, 
# 		since removing orthologs with low values can improve tree certainty
# 30 OCT 2017: requires job type in command, to avoid confusion when running with wrong job type set in script
# 
# 
# 
# 
# 
#######################################################################


$tree1 = $ARGV[0];
$which_job = $ARGV[1];

unless($tree1 =~ /[\w\d]/ && $which_job =~ /\-bootstrap|\-var/)
	{
	die "\n\tCOMMAND ERROR. give file name and job type. quitting.\n"	
	};
if($which_job =~ /boot/)
	{
	$calculate_average_bootstrap_support = 1;
	print "this script calculates either average_bootstrap or branch length variance, you have opted for FORMER\n";
	}else{
	print "this script calculates either average_bootstrap or branch length variance, you have opted for LATTER\n";
	};


$bootstrap_cutoff = 0.70;

$output = "all_branch_results";# all_branch_results_raxml prev: all_branch_results










# writes file called varianceR

##################
write_R_script();#
##################



	#######################
$newick = read_tree($tree1);#
	#######################


#print "tree1 start:$newick\n";


	########################################
$newick1 = standarise_branchlengths($newick);#
	########################################



if ($calculate_average_bootstrap_support == 1)
	{

	my $sum_boot = 0;my $count_boot = 0;
#	print "\ntestA newick1:$newick1\n";
	while($newick1 =~ s/\)([\d\.]+)/)/)
		{
		my $boot = $1;$sum_boot += $boot;$count_boot++;
		#print "$boot ";
		};
#	print "\n\ntestB newick1:$newick1\n";
	my $av_boot = $sum_boot / $count_boot;	
	print "
input file:$tree1
sum_boot:$sum_boot count_boot:$count_boot
average bootstrap support:$av_boot\n\n";

	if($av_boot >= $bootstrap_cutoff)
		{
		my $fasta_file = $tree1;
		unless($fasta_file =~ s/(insectaNUCL.+.clo_pruned.fa).orth_ft_pruned/$1/){die "\nerror 89.\n"};
		system("cp $fasta_file $fasta_file.RM5");
		print "av boot ($av_boot) >= your cutoff $bootstrap_cutoff, KEEPING ortholog\n";
		}else{
		print "av boot ($av_boot) < your cutoff $bootstrap_cutoff, discarding ortholog\n";
		};


	print "\n\$calculate_average_bootstrap_support == $calculate_average_bootstrap_support\n" , 
		"not calculating branchlength var\nhave written new file $fasta_file.RM5\n ... quitting\n";
	exit;	

	};# if $calculate_average_bootstrap_support == 1






my @blengths = ();
my $sum_branchlength =0;
my $count_bls = 0;
open(OUT, ">table_for_R") || die "\nerror 3\n";

# print "newickA:($newick1)\n";

# presumably this was developed on fasttree tree, i just checked and it works perfect on raxml trees also

while($newick1 =~ /(\:)(\d+\.\d+)([\(\)\,])/)
	{
	my $len = $2;my $following_char = $3;
	print OUT "$len\n";# print "$len\n";
	$sum_branchlength += $len;
	$count_bls++;
	push @blengths , $len;
	#$len *= 0.1;
	$newick1 =~ s/(\:)(\d+\.\d+)([\(\)\,])/$3/; 
	}

close OUT;

# print "newickB:($newick1)\n";



$cooamdn = "R --vanilla < varianceR";
print "running command:$cooamdn\n";



system ($cooamdn);


open(RES , "R_Variance_OUT");# || die "\nprint 58\n";
my $sum = "NA";my $vari = "NA";
while (my $line = <RES>)
	{
	if($line =~ /^([\.\d]+)\s([\.\d]+)/)
		{
		$sum = $1; $vari = $2;
		};
	}
close RES;

open(RES2 , ">>$output") || die "";
print RES2 "$tree1\t$sum\t$vari\n";
close RES2;

print  "results:\n$tree1 $sum $vari
output file is named $output
\n";


print "\nFIN.\n";
exit;

##########################################################################################################


sub standarise_branchlengths
{

my $current_tree = shift;

# ITS1YS33:-0.00000002)	replace negative branchlength with zero
if($current_tree =~ s/\:\-\d+\.\d+(\D)/\:0\.0$1/g)
	{
	print "Negative branchlengths found and replaced\n";
	}else{
	print "passed check 1\n";
	};


# ITS1YS33:3.892891902e-06, replace scientific notation
if($current_tree =~ s/\:(\d\.\d+)[eE](\-\d+)/"\:" .  ($1 * ( 10 ** $2 ))/ge)
	{
	print "branchlengths in scientific notation found, replaced\n";
	}else{
	print "passed check 2\n";
	};




return($current_tree);

}





##########################################################################################################


sub read_tree
{
my $file = shift;

open(IN, $file) || die "error 206, cant open file:($file)\n";
my $t = "";my $found=0;
while(my $line = <IN>)
	{
	if($line =~ /^(.+\(.+)$/)
		{
		$t = $1;$found=1;
		}
	}
close(IN);
unless($found == 1){die "cant find tree\n";}
print "got tree from file $file\n";

return($t);


};




##########################################################################################################



sub write_R_script
{

open(OUT_R, ">varianceR") || die "\nerror cant open out:varianceR\n";

print OUT_R "

t1<-read.table(\"table_for_R\", header=F);
sum1<-sum(t1)
varience1<-var(t1)
print (c(\"sum \" , sum1 , \"varience \" , varience1))
write(c(sum1 , varience1), file = \"R_Variance_OUT\",
      ncolumns = 2, append = FALSE, sep = \" \")

";

close OUT_R;

};



##########################################################################################################




