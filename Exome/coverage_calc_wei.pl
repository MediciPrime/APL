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
	$cmd = q(awk '{if($NF > 0) print}');
	system "bedtools coverage -split -abam $1.sorted.rmdup.bam -b /data/679655_LH_Leidos/Exome/illumina.bed | $cmd > $1.Illumina.exon.coverage.bed";
	system "bedtools coverage -split -abam $1.sorted.rmdup.bam -b /data/reference/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf | $cmd > $1.UCSC.exon.coverage.gtf";
	chdir "../..";
	print "***coverage_cal finished:	$1\n";
}