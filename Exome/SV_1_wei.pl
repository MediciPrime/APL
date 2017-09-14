#!/usr/bin/perl -w

# by Wenhu Guo (wenhu_guo@acgtinc.com)
#

use Cwd;

foreach(glob("*_D[MN]*/analysis/*.sorted.rmdup.bam")){
	chomp;
	$_ =~ m/(\d+_D[MN].*)\/analysis\/(\d+_D[MN].*)\.sorted.rmdup.bam/;
	if ("$1" ne "$2") {
	   die "Folder name and file name don't match! Fatal error.\n";
	}
	chdir "$1/analysis";
	system "mkdir SV";
	system "/data/bin/krishna/Leidos/samtools-0.1.19/samtools view $2.sorted.rmdup.bam | tail -n+1000000 | /data/bin/krishna/Leidos/lumpy-sv/scripts/pairend_distro.py -r 150 -X 4 -N 10000 -o SV/$2.hist 1> SV/$2.stats";
	chdir "../..";
	print "***SV1 finished:	$2\n";
}