#!/usr/bin/perl -w

# by Wei Zhu (wei_zhu@acgtinc.com)
#

foreach(glob("*_D[MN]*/analysis/*.exome.vcf")){
	chomp;
	$_ =~ m/(\d+_D[MN].*)\/analysis\/(\d+_D[MN].*)\.exome\.vcf/;
	if ("$1" ne "$2") {
	   die "Folder name and file name don't match! Fatal error.\n";
	}
	system "python2.7	/data/bin/krishna/Leidos/annotate.py	$_	wei_zhu@acgtinc.com";
	print "Done: $_\n";
}