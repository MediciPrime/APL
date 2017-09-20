#!/usr/bin/perl -w

#
# by  Wei Zhu (wei_zhu@acgtinc.com)
# 
# Created date: 04/11/2016 

use strict;

# Provide project name and input data
die "usage:$0 <project_name>\t<reference_sequence>\t<num_of_threads>\n" unless $#ARGV == 2;
my $project_name = shift;
my $ref_seq = shift;
my $num_of_threads = shift;

system "samtools view -h -b -F 4 $project_name.sam > $project_name.bam";

system "samtools sort -@ $num_of_threads $project_name.bam $project_name.sorted";
system "rm -f $project_name.bam";
system "samtools index $project_name.sorted.bam";
system "samtools depth $project_name.sorted.bam > $project_name.depth";

system "samtools faidx $ref_seq";
# reweigh BAQ scores to reduce false positive
system "samtools calmd -Abr $project_name.sorted.bam $ref_seq > $project_name.tmp";
system "mv $project_name.tmp $project_name.sorted.bam";
system "samtools index $project_name.sorted.bam";
system "samtools mpileup -d 1000000 -f $ref_seq $project_name.sorted.bam > $project_name.mpileup";
system "cut -f1-4 $project_name.mpileup > $project_name.real.depth";