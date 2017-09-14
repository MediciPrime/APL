#!/usr/bin/perl -w

# by Wenhu Guo (wenhu_guo@acgtinc.com)
#

use Cwd;

foreach(glob("*_D[MN]*/fastq/*.R1_val_1.fq.gz")){
	chomp;
	$_ =~ m/(\d+_D[MN].*)\/fastq\/(\d+_D[MN].*)\.R1_val_1\.fq\.gz/;
	if ("$1" ne "$2") {
	   die "Folder name and file name don't match! Fatal error.\n";
	}
	$read1 = $2.".R1_val_1.fq.gz";
	$read2 = $2.".R2_val_2.fq.gz";
	chdir "$1/analysis";
	#print(cwd);
	system "/data/bin/krishna/Leidos/bwa-0.7.8/bwa mem -M -T 45 -t 56 /data/reference/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa ../fastq/$read1 ../fastq/$read2 > $2.sam";
	system "/data/bin/krishna/Leidos/samtools-0.1.19/samtools view -@ 56 -h -bS $2.sam > $2.bam";
	system "/data/bin/krishna/Leidos/samtools-0.1.19/samtools sort -@ 56 $2.bam $2.sorted";
	system "/data/bin/krishna/Leidos/samtools-0.1.19/samtools index $2.sorted.bam";
	system "rm $2.sam";
	system "rm $2.bam";
	system "java -Xmx10g -jar /data/bin/krishna/Leidos/picard-tools-1.126/picard.jar MarkDuplicates I=$2.sorted.bam REMOVE_DUPLICATES=TRUE O=$2.sorted.rmdup.bam M=$2.rmdup.metrics";
	chdir "../..";
	#print(cwd);
}