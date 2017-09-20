#!/usr/bin/perl -w

#
# by  Wei Zhu (wei_zhu@acgtinc.com)
# 

use strict;

# Provide project name and input data
die "usage:$0 <depth_file>\t<contig_prefix: eg. scaffold>\t<num_of_contigs: defualt 50>\n" unless $#ARGV == 1 or $#ARGV == 2;
my $depth_file = shift;
my $contig_prefix = shift;
my $num_contigs = 50;
$num_contigs = shift if ( $#ARGV == 2 );

foreach my $i (1..$num_contigs) {
	my $contig = $contig_prefix.$i;
	#system "grep $contig $depth_file -n -m 1 -o"; 
	system "grep $contig $depth_file -n -m 1 | cut -f1 -d:"; 
}
