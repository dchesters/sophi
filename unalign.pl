


# 
# required if you want to make a blast database out of aligned seqs,
# gaps need removing
# 
# 
# 
# 
# 
# CHANGE LOG 
# 
# 2020-01-15: prints entries count, to avoid confusion
# 2020-02-10: doesnt print line if contains nothing after gap removal. (this happens a lot due to some linux idiosyncrasy somewhere)
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
# 
# 
# 
###################################################################################################################

$in = $ARGV[0];

open(IN, $in) || die "\n";
open(OUT, ">$in.unaligned") || die "\n";


while(my $line = <IN>)
	{
	if($line =~ /^>/)
		{
		print OUT $line;$fasta_entries++;
		}else{
		if($line =~ s/\-//g){$gaps_removed_from_seqs++};

		if($line =~ /./){print OUT $line};
		
		}

	}


close IN;
close OUT;


print "
file name:$in
contains entries count:$fasta_entries
total gaps removed from file:$gaps_removed_from_seqs
";


