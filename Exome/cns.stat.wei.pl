#!/usr/bin/perl -w

# by Wei Zhu (wei_zhu@acgtinc.com)
#

foreach(glob("*_D[MN]*/analysis/*.exome.cns.fasta")){
	chomp;
	$_ =~ m/(\d+_D[MN].*)\/analysis\/(\d+_D[MN].*)\.exome\.cns\.fasta/;
	if ("$1" ne "$2") {
	   die "Folder name and file name don't match! Fatal error.\n";
	}
	system "assemblathon_stats.pl  $_ | head -n 23 | sed 's/scaffold/exon/g' > $2/analysis/$2.exome.cns.stats";
	print "Done: $_ --> $2.exome.cns.stats\n";
}