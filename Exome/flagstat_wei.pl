#!/usr/bin/perl -w

# by Wenhu Guo (wenhu_guo@acgtinc.com)
#

use Cwd;

foreach(glob("*_D[MN]*/analysis/*.sorted.bam")){
	chomp;
	$_ =~ m/(\d+_D[MN].*)\/analysis\/(\d+_D[MN].*)\.sorted\.bam/;
	if ("$1" ne "$2") {
	   die "Folder name and file name don't match! Fatal error.\n";
	}
	chdir "$1/analysis";
	system "samtools flagstat $2.sorted.bam > $2.samtools.flagstat.out";
	chdir "../..";
	print "***flagstat finished:	$2\n";
}