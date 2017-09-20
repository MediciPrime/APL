#!/usr/bin/perl -w

#
# extract sequence reads that contain motifs from fastq file
# by Wenhu Guo (wguo@unl.edu)
# Last update: 06/05/2013
#

use strict;

die "usage: $0 <Solexa-PE-fastq-File>\t<OUT FILE>\t<motif SEQ>\n" unless $#ARGV == 2;
my $SolexaFile = shift;
my $outfile =shift;
my $motif = shift;
my $motif_rc = reverse $motif;
$motif_rc =~ tr/ACGTURYMKWSBDHVNacgturymkwsbdhvn/TGCAAYRKMWSVHDBNtgcaayrkmwsvhdbn/;

open (IN, "$SolexaFile") or die "can't open $SolexaFile\n";
open (OUT, ">$outfile") or die "can't open file\n";
my $ReadID = "";
my $ReadSeq = "";
my $Read_qual_flag = "";
my $ReadQuality = "";
my $NeedSeq = 0;
my $NeedQual = 0;
while (<IN>) {
    chomp;
    $ReadID = $_;
	$ReadSeq = <IN>;
	$Read_qual_flag = <IN>;
	$ReadQuality = <IN>;
	if ($ReadSeq =~ m/($motif)+/gi) {
	   print OUT ">"."$ReadID\n"."$ReadSeq\n";
    } elsif ($ReadSeq =~ m/($motif_rc)+/gi) {
	   my $ReadSeqRC = reverse $ReadSeq;
	   $ReadSeqRC =~ tr/ACGTURYMKWSBDHVNacgturymkwsbdhvn/TGCAAYRKMWSVHDBNtgcaayrkmwsvhdbn/;
	   print OUT ">"."$ReadID\n"."$ReadSeqRC\n";
	}
}
close OUT;
close IN;
