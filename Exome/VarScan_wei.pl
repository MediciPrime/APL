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
	system "/data/bin/krishna/Leidos/samtools-0.1.19/samtools index $2.sorted.rmdup.bam";
	system "/data/bin/krishna/Leidos/samtools-0.1.19/samtools mpileup -l /data/679655_LH_Leidos/Exome/illumina.bed -f /data/reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa $2.sorted.rmdup.bam | VarScan mpileup2cns --variants 1 --output-vcf 1 --min-avg-qual 25 > $2.exome.vcf";
	chdir "../..";
	print "****VarScan finished:	$2\n";
}