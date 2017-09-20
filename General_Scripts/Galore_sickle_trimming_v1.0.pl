#!/usr/bin/perl -w
#
# Common usage: Galore_sickle_trimming.pl project_name sample.R1.fq sample.R2.fq 56 1 1

use strict;

# Provide project name and input data
die "usage:$0 <project_name>\t<read1_forward>\t<read2_reverse>\t<num_of_threads>\t<trim_galore_trimming: 0 or 1>\t<sickle_trimming: 0 or 1>" unless $#ARGV == 5;
my $project_name = shift;
my $read1 = shift;
my $read2 = shift;
my $num_of_threads = shift;
my $trim_galore_trimming = shift;
my $sickle_trimming = shift;
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
   system "trim_galore -q 30 --length 50 --no_report_file --paired $read1 $read2 2>$project_name.trim_galore.err";
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