#!/usr/bin/perl -w

# by Wenhu Guo (wenhu_guo@acgtinc.com)
#

foreach(glob("*.fq")){
	chomp;
	$in = $_;
	$in =~ s/\.fq//;
	system "cutadapt -g GTTCAGAGTTCTACAGTCCGACGATC -a TGGAATTCTCGGGTGCCAAGG -q 30 -m 17 -M 26 -o $in.trimmed.fq $_ 1>$in.cutadapt.log";
	$cmd_tmp = q(awk '{if(NR%4==2) print length($1)}');
	system "cat $in.trimmed.fq | $cmd_tmp > $in.trimmed.readsLength.txt\n";
    system "./textHistogram -maxBinCount=75 $in.trimmed.readsLength.txt > $in.trimmed.readsLength.freq\n";
	system "./groupReads.pl input=$in.trimmed.fq output=$in.trimmed.rc";
	system "bowtie2 -x /data/reference/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome -D 20 -R 3 -N 1 -L 15 -i S,1,0.50 -U $in.trimmed.fq -p 45 -S $in.sam --al $in.trimmed.mapped.fq --un $in.trimmed.unmapped.fq 2>$in.bowtie2.err";
    system "samtools view -bS -@ 45 -o $in.bam $in.sam";
    system "samtools sort -@ 45 $in.bam $in.sorted";
	system "samtools index $in.sorted.bam";
}