#!/usr/bin/perl -w

#
# by  Wenhu Guo (wenhu_guo@acgtinc.com)
# 
# Created date: 09/30/2011 
#
# An integrated pipeline for de novo genome sequencing project (for plasmid to bacteria-sized genomes)
# 1) trim raw reads (trim_galore) 
# 2) de novo assembly (SPAdes,Velvet) 
# 
#
# Common usage: de_novo_sequencing_pipeline_v1.0.pl project_name sample.R1.fastq sample.R2.fastq 31 91 10 30 150 30 10 600 0 0
# Common usage: de_novo_sequencing_pipeline_v1.4.pl project_name sample.R1.fastq sample.R2.fastq 71 131 10 300 600 100 10 550 1 1 0 1 1 0.001
# Common usage: de_novo_sequencing_pipeline_v1.5.pl project_name sample.R1.fastq sample.R2.fastq 1 1 0 1 
# Version update:
# 1.01, 09/23/2014, fixed the output of comparison_of_assemblies.txt
# 1.10, 11/10/2014, added the command line of spade assembly, print the progress to STDOUT 
# 1.11, 1/12/2015, changed the output names
# 1.20, 1/15/2015, changed the parameters for spade assembly
# 1.21, 1/30/2015, no need to combine R1 and R2 reads for velvet
# 1.22, 2/16/2015, provided options to choose an assembler
# 1.30, 2/16/2015, provided an option to use a subset of PE reads
# 1.31, 4/7/2015, updated output folder names
# 1.40, 11/16/2015, added an option of sickle trimming <sickle_trimming: 0 or 1>; changed <trimmed_reads: 0  or 1> to <trim_galore_trimming: 0 or 1>; changed "read[1|2]" to "read[1|2]_trimmed" in all steps following trimming -Wei Zhu
# 1.50, 09/06/2016, velvet assembly was no longer needed by deafult

use strict;


die "usage:$0 <project_name>\t<reads_R1>\t<reads_R2>\t<trim_galore_trimming: 0 or 1>\t<sickle_trimming: 0 or 1>\t<mate_pair: 0 or 1>\t<Fraction_reads_to_use: 0 to 1>\n" unless $#ARGV == 6;

my $project_name = shift;
my $read1 = shift;
my $read2 = shift;
my $trim_galore_trimming = shift;
my $sickle_trimming = shift;
my $mate_pair = shift;
my $FracToKeep = shift;
my $read1_trimmed;
my $read2_trimmed;
my $single_reads;

# if reads were already trimmed
unless ($trim_galore_trimming || $sickle_trimming) {
       $read1_trimmed = $read1;
	   $read2_trimmed = $read2;
}

# trim reads with trim_galore
if ($trim_galore_trimming) {
   print "Start trim_galore trimming: q 30, length 50...\n";
   system "trim_galore -q 30 --length 50 --paired $read1 $read2 2>$project_name.trim_galore.err";
   if ($read1 =~ /\.fastq$/ && $read2 =~ /\.fastq$/){
       $read1 =~ s/\.fastq$/_val_1.fq/;
	   $read2 =~ s/\.fastq$/_val_2.fq/;
   } elsif ($read1 =~ /\.fq$/ && $read2 =~ /\.fq$/){
       $read1 =~ s/\.fq$/_val_1.fq/;
	   $read2 =~ s/\.fq$/_val_2.fq/;
   } else{
       $read1 =~ s/\.fastq$/_val_1.fq/;
	   $read2 =~ s/\.fastq$/_val_2.fq/;
   }
   $read1_trimmed = $read1;
   $read2_trimmed = $read2;
}

# trim reads with sickle (necessary if raw reads have low quality scores)
if ($sickle_trimming) {
   if ($trim_galore_trimming) {
	  print "Start sickle_trimming: q 30, length 50...\n";	
      my $prefix1 = $read1_trimmed;
	  my $prefix2 = $read2_trimmed;
	  $prefix1 =~ s/_val_1\.fq//;
	  $prefix2 =~ s/_val_2\.fq//;
	  $read1_trimmed = $prefix1."_val_1.sickle.fq";
	  $read2_trimmed = $prefix2."_val_2.sickle.fq";
	  $single_reads = $prefix1.".sickle.single.fq";
      system "sickle pe -q 30 -l 50 -t sanger -f $read1 -r $read2 -o $read1_trimmed -p $read2_trimmed -s $single_reads 1>$project_name.sickle.log";
	}
}


if ($FracToKeep <= 0 or $FracToKeep > 1) {
   die "Fraction of reads to use must be great than 0 and less or equal to 1\n";
} elsif ($FracToKeep < 1) {
   system "PullOutPEReadSubset_v1.0.pl $read1_trimmed $read2_trimmed $project_name $FracToKeep";
   $read1_trimmed = "$project_name.R1.$FracToKeep.fq";
   $read2_trimmed = "$project_name.R2.$FracToKeep.fq";
}


print "Start Spades assembly...\n";
    if ($mate_pair) {
        system "spades.py --hqmp1-1 $read1_trimmed --hqmp1-2 $read2_trimmed --hqmp1-rf --careful -o $project_name.$FracToKeep.spades_assembly";
    } else {
        system "spades.py --pe1-1 $read1_trimmed --pe1-2 $read2_trimmed --careful -o $project_name.$FracToKeep.spades_assembly";
    }
