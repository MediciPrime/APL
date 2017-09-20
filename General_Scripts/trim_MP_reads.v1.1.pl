#!/usr/bin/perl -w    

# by Wenhu Guo (wenhu_guo@acgtinc.com)
# This script was used to trim mate-pair reads
# Please refer to this technical note, http://www.illumina.com/documents/products/technotes/technote_nextera_matepair_data_processing.pdf

use strict;

# Provide project name and input data
die "usage:$0 <project_name>\t<read1_forward>\t<read2_reverse>\n" unless $#ARGV == 2;

my $project_name = shift;
my $read1 = shift;
my $read2 = shift;

# first trim Nextera Transposase adapter 
system "cutadapt -b CTGTCTCTTATA -B CTGTCTCTTATA -m 35 -q 20 -o $project_name.R1.trim1.fastq -p $project_name.R2.trim1.fastq $read1 $read2 1>$project_name.trim1.log";

# second trim external adapter 
system "cutadapt -g GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -G GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -m 35 -q 20 -o $project_name.R1.trim2.fastq -p $project_name.R2.trim2.fastq $project_name.R1.trim1.fastq $project_name.R2.trim1.fastq 1>$project_name.trim2.log";

# then trim junction adapter
system "cutadapt -b CTGTCTCTTATACACATCT -b AGATGTGTATAAGAGACAG -B CTGTCTCTTATACACATCT -B AGATGTGTATAAGAGACAG -m 35 -q 20 -o $project_name.R1.trimmed.fastq -p $project_name.R2.trimmed.fastq $project_name.R1.trim2.fastq $project_name.R2.trim2.fastq 1>$project_name.trim3.log";
