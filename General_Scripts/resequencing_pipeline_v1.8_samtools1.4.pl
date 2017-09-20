#!/usr/bin/perl -w

#
# by  Wenhu Guo & Wei Zhu (wei_zhu@acgtinc.com)
# 
# Created date: 02/24/2012 
# Last update: 09/22/2014, version 1.00
#              10/07/2014, added the option of sickle trimming
#			   10/07/2014, Recalc BAQ (/data/bin/samtools-1.4/samtools calmd -Abr) scores using /data/bin/samtools-1.4/samtools to reduce false positives in BAM file -- Krishna Motheramgari
#              10/13/2014, fixed the bug of filter out non-aligning reads, also use /data/bin/samtools-1.4/samtools to call variants.
#              11/12,2014, changed the parameters for VarScan
#              11/16/2014, updated the maximum mpileup depth from 8000 to 1000000
#               1/21/2015, changed the name of a temporary file
#				1/28/2015, provided an option of a quick bowite alignment (no variant calling) 
#               2/10/2015, updated the command lines for variant calling and consensus calling 
#	            3/20/2015, new index file is generated after calmd
#               5/1/2015, updated the parameter for VarScan (from --strand-filter 0 to default)
#               5/4/2015, added an option for using bwa for stringent mapping (MAPQ 30)
#		        9/22/2016, added freebayes for variant calling and an option to use local mapping with bowtie2, version 1.6 -- Wei Zhu
#				04/18/2017, update samtools and bcftools to v1.4; added --ploidy to freebayes and --aa (absolute all positions) & -L (max per-sample depth for INDEL calling) to samtools. v1.7 -- Wei Zhu
#				05/10/2017, update trim_galore and sickle to handel gzip'ed fastq input and output. v1.8 -- Wei Zhu 
# An integrated pipeline for genome resequencing project
# 1) trim raw reads (trim_galore and/or sickle) 
# 2) map trimmed reads to reference sequence (bowtie2 or BWA)
# 3) generate sorted BAM file (/data/bin/samtools-1.4/samtools) 
# 4) call variants against reference sequence (freebayes, VarScan and /data/bin/samtools-1.4/samtools)
# 5) generate consensus sequence (/data/bin/samtools-1.4/samtools)
#
# Common usage: resequencing_pipeline.pl project_name reference.fas sample.R1.fq sample.R2.fq 0 1000 16 1 0 0 0 1

use strict;

# Provide project name and input data
die "usage:$0 <project_name>\t<reference_sequence>\t<read1_forward>\t<read2_reverse>\t<insert_size_min>\t<insert_size_max>\t<num_of_threads>\t<trim_galore_trimming: 0 or 1>\t<sickle_trimming: 0 or 1>\t<stringent mapping: 0 or 1>\t<local mapping: 0 or 1>\t<variant_calling: 0 or 1>\n" unless $#ARGV == 11;
my $project_name = shift;
my $ref_seq = shift;
my $read1 = shift;
my $read2 = shift;
my $insert_min = shift;
my $insert_max = shift;
my $num_of_threads = shift;
my $trim_galore_trimming = shift;
my $sickle_trimming = shift;
my $stringent_mapping =shift;
my $local_mapping =shift;
my $variant_calling = shift;
my $ConsensusSeq = "$project_name".".consensus.fas";
my $read1_trimmed;
my $read2_trimmed;
my $single_reads;
my $depth = 10000000;
my $ploidy = 1;

# if reads were already trimmed
unless ($trim_galore_trimming || $sickle_trimming) {
       $read1_trimmed = $read1;
	   $read2_trimmed = $read2;
}

# trim reads with trim_galore
if ($trim_galore_trimming) {
   system "trim_galore -q 30 --length 50 --no_report_file --paired $read1 $read2 2>$project_name.trim_galore.err";
   if ($read1 =~ /\.fastq.gz$/ && $read2 =~ /\.fastq.gz$/){
       $read1 =~ s/\.fastq.gz$/_val_1.fq.gz/;
	   $read2 =~ s/\.fastq.gz$/_val_2.fq.gz/;
   } elsif ($read1 =~ /\.fq.gz$/ && $read2 =~ /\.fq.gz$/){
       $read1 =~ s/\.fq.gz$/_val_1.fq.gz/;
	   $read2 =~ s/\.fq.gz$/_val_2.fq.gz/;
   } else{
       $read1 =~ s/\.fastq.gz$/_val_1.fq.gz/;
	   $read2 =~ s/\.fastq.gz$/_val_2.fq.gz/;
   }
   $read1_trimmed = $read1;
   $read2_trimmed = $read2;
}

# trim reads with sickle (necessary if raw reads have low quality scores)
if ($sickle_trimming) {
   if ($trim_galore_trimming) {
      my $prefix1 = $read1_trimmed;
	  my $prefix2 = $read2_trimmed;
	  $prefix1 =~ s/_val_1\.fq.gz//;
	  $prefix2 =~ s/_val_2\.fq.gz//;
	  $read1_trimmed = $prefix1."_val_1.sickle.fq.gz";
	  $read2_trimmed = $prefix2."_val_2.sickle.fq.gz";
	  $single_reads = $prefix1.".sickle.single.fq.gz";
      system "sickle pe --gzip-output -q 30 -l 50 -t sanger -f $read1 -r $read2 -o $read1_trimmed -p $read2_trimmed -s $single_reads 1>$project_name.sickle.log";
	}
}
    
# mapping trimmed reads to reference sequence
if ($stringent_mapping) {
    system "bwa index -p $project_name $ref_seq";
	system "bwa mem -t $num_of_threads -T 30 $project_name $read1_trimmed $read2_trimmed > $project_name.sam 2>$project_name.bwa.err";
	# filter out non-aligning reads
	system "/data/bin/samtools-1.4/samtools view -h -b -F 4 -@ $num_of_threads $project_name.sam > $project_name.bam";	
} else {
	system "bowtie2-build $ref_seq $project_name";
	if ($local_mapping) {
		system "bowtie2 --local -x $project_name -1 $read1_trimmed -2 $read2_trimmed -I $insert_min -X $insert_max --fr -p $num_of_threads -S $project_name.sam 2>$project_name.bowtie2.err";
	}
	else {
		system "bowtie2 -x $project_name -1 $read1_trimmed -2 $read2_trimmed -I $insert_min -X $insert_max --fr -p $num_of_threads -S $project_name.sam 2>$project_name.bowtie2.err";
	}
	# filter out non-aligning reads
	system "/data/bin/samtools-1.4/samtools view -h -b -F 4 -@ $num_of_threads $project_name.sam > $project_name.bam";
}

system "/data/bin/samtools-1.4/samtools sort -@ $num_of_threads $project_name.bam -o $project_name.sorted.bam";
system "rm -f $project_name.bam";
system "/data/bin/samtools-1.4/samtools index -@ $num_of_threads $project_name.sorted.bam";
system "/data/bin/samtools-1.4/samtools depth -aa $project_name.sorted.bam > $project_name.depth";

# skip variant calling for a quick bowtie/bwa alignment
if ($variant_calling) {
    system "/data/bin/samtools-1.4/samtools faidx $ref_seq";
    # reweigh BAQ scores to reduce false positive
    system "/data/bin/samtools-1.4/samtools calmd -Abr -@ $num_of_threads $project_name.sorted.bam $ref_seq > $project_name.tmp";
    system "mv $project_name.tmp $project_name.sorted.bam";
	system "/data/bin/samtools-1.4/samtools index -@ $num_of_threads $project_name.sorted.bam";
	system "/data/bin/samtools-1.4/samtools mpileup -d $depth -aa -f $ref_seq $project_name.sorted.bam > $project_name.mpileup";
	system "cut -f1-4 $project_name.mpileup > $project_name.real.depth";
    
	# call variants (SNPs and INDELs) with VarScan (non-default parameters)
    system "VarScan mpileup2cns $project_name.mpileup --variants 1 --min-coverage 10 --min-reads2 3 --min-avg-qual 30 --min-var-freq 0.05 > $project_name.variants.freq0.05.out";
    # reformat the output of variant information
    system "VarParse.pl <$project_name.variants.freq0.05.out > $project_name.variants.freq0.05.txt";
    
	# call variants with samtools & bcftools multiallelic calling model 
    system "/data/bin/samtools-1.4/samtools mpileup -d $depth -L $depth -ugf $ref_seq $project_name.sorted.bam | /data/bin/bcftools-1.4/bcftools call --ploidy $ploidy -vm -O v -o $project_name.samtools.vcf";
	
	# call variants with samtools & bcftools consensus calling model 
	system "/data/bin/samtools-1.4/samtools mpileup -d $depth -L $depth -uf $ref_seq $project_name.sorted.bam | /data/bin/bcftools-1.4/bcftools call --ploidy $ploidy -vc -o $project_name.samtools.consensus.vcf";
    
	# generate consensus sequence with /data/bin/samtools-1.4/samtools
    system "/data/bin/samtools-1.4/samtools mpileup -d $depth -L $depth -uf $ref_seq $project_name.sorted.bam | /data/bin/bcftools-1.4/bcftools call --ploidy $ploidy -c | vcfutils_indel_ACGT_v1.0.pl vcf2fq -d 100 -D $depth > $ConsensusSeq";
	
	# call variants with freebayes
    system "freebayes -f $ref_seq --ploidy $ploidy -F 0.05 --min-coverage 10 $project_name.sorted.bam > $project_name.freebayes.vcf";
    system "vcffilter -f 'QUAL > 20' $project_name.freebayes.vcf > $project_name.freebayes.Q20.vcf";
	

}
