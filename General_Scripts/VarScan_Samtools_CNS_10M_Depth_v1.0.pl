#!/usr/bin/perl -w
# Call variants using both VarScan and samtools from sorted bam file with depth up to 10 million. Also generate consensuse sequences.  

die "usage:$0 <project_name>\t<reference_sequence>\n" unless $#ARGV == 1;
my $project_name = shift;
my $ref_seq = shift;
my $ConsensusSeq = "$project_name".".consensus.fas";

system "samtools mpileup -d 10000000 -f $ref_seq $project_name.sorted.bam > $project_name.mpileup";
system "cut -f1-4 $project_name.mpileup > $project_name.real.depth";
# call variants (SNPs and INDELs) with VarScan (non-default parameters)
system "VarScan mpileup2cns $project_name.mpileup --variants 1 --min-coverage 10 --min-reads2 3 --min-avg-qual 30 --min-var-freq 0.05 > $project_name.variants.freq0.05.out";
# reformat the output of variant information
system "VarParse.pl <$project_name.variants.freq0.05.out > $project_name.variants.freq0.05.txt";
# call variants with SAMtools
system "samtools mpileup -d 10000000 -ugf $ref_seq $project_name.sorted.bam | bcftools call -vm -O z -o $project_name.SAMtools.vcf.gz";
# generate consensus sequence with SAMtools
system "samtools mpileup -d 10000000 -uf $ref_seq $project_name.sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > $ConsensusSeq";
