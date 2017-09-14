#!/usr/bin/perl -w

#
# by  Wei Zhu (wei_zhu@acgtinc.com)
#

use strict;

die "usage:$0 <read1_forward>n" unless $#ARGV == 0;

my $read1 = shift;
my $read2 = $read1;

if ($read1 =~ m{R1\.fq\.gz$}) {
	$read2 =~ s{R1}{R2}; 
}
else {
	print "\nERROR!!! Input must be sample_R1.fq.gz!!!\n\n";
	exit;
}

print "Start	trim_galore $read1 <------> $read2\n";
system "trim_galore -q 30 --length 50 --suppress_warn --paired $read1 $read2";
print "Done	trim_galore $read1 <------> $read2\n";
