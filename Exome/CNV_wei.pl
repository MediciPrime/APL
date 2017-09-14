#!/usr/bin/perl -w

# by Wenhu Guo (wenhu_guo@acgtinc.com)
#

use Cwd;

foreach(glob("*_DM*/analysis/*.sorted.rmdup.bam")){
	chomp;
	$_ =~ m/(\d+_DM.*)\/analysis\/(\d+_DM.*)\.sorted.rmdup.bam/;
	if ("$1" ne "$2") {
	   die "Folder name and file name don't match! Fatal error.\n";
	}
	chdir "$1/analysis";
	$prefix = $1;
	$DM_sample_name = $2;
	$prefix =~ s/\_DM//;
	$DN_sample_name = "$prefix"."_DN";
	system "/data/bin/krishna/Leidos/CONTRA.v2.0.6/contra.py -f /data/reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -t /data/679655_LH_Leidos/Exome/illumina.bed -c /data/679655_LH_Leidos/Exome/Batch7/cat_set/p1/$DN_sample_name/analysis/$DN_sample_name.sorted.rmdup.bam -s $DM_sample_name.sorted.rmdup.bam -o CNV";
	chdir "../..";
	print "Done: $_\n";
}