#!/usr/bin/perl -w

#
# extract sequence reads that contain motifs from fasta file
# by Wenhu Guo (wguo@unl.edu)
# Last update: 06/05/2013
#

use strict;

die "usage: $0 <FASTA-File>\t<OUT FILE>\t<motif SEQ>\n" unless $#ARGV == 2;
my $FastaFile = shift;
my $outfile = shift;
my $motif = shift;

open (IN, "$FastaFile") or die "can't open $FastaFile\n";
open (OUT, ">$outfile") or die "can't open file $outfile\n";
my $ReadID = "";
my $ReadSeq = "";
my $motif_rc = reverse $motif;
$motif_rc =~ tr/ACGTURYMKWSBDHVNacgturymkwsbdhvn/TGCAAYRKMWSVHDBNtgcaayrkmwsvhdbn/;

while (<IN>) {
    chomp;
	$ReadID = $_;
    $ReadSeq = <IN>;
	if ($ReadSeq =~ m/($motif)+/gi) {
	   print OUT "$ReadID\n"."$ReadSeq\n";
    } elsif ($ReadSeq =~ m/($motif_rc)+/gi) {
	   my $ReadSeqRC = reverse $ReadSeq;
	   $ReadSeqRC =~ tr/ACGTURYMKWSBDHVNacgturymkwsbdhvn/TGCAAYRKMWSVHDBNtgcaayrkmwsvhdbn/;
	   print OUT "$ReadID\n"."$ReadSeqRC\n";
	}
}
close OUT;
close IN;
