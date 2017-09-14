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
	system "/data/bin/krishna/Leidos/samtools-0.1.19/samtools mpileup -l /data/679655_LH_Leidos/Exome/illumina.bed -uf /data/reference/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa $2.sorted.rmdup.bam | /data/bin/krishna/Leidos/samtools-0.1.19/bcftools/bcftools view -cg -| /data/bin/krishna/Leidos/samtools-0.1.19/bcftools/vcfutils.pl vcf2fq > $1.cns.fasta";
	system "bedtools getfasta -split -fi $1.cns.fasta -bed $1.Illumina.exon.coverage.bed -fo $1.exome.cns.fasta";
	system "/data/bin/krishna/Leidos/trans < $1.exome.cns.fasta > $1.exome.cns.aa";
	chdir "../..";
	print "***consensus_calling finished:	$1\n";
}