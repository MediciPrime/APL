#!/usr/bin/perl -w

#
# extract sequence reads that contain both 3'end and 5'end motifs
# by Wenhu Guo (wenhu_guo@acgtinc.com)
# 02/09/2015
#

use strict;

die "usage: $0 <Fastq-File>\t<OUT-FILE>\t<motif1-SEQ>\t<motif2-SEQ>\n" unless $#ARGV == 3;
my $FastqFile = shift;
my $outfile =shift;
my $motif1 = shift;
my $motif2 = shift;
my $motif1_rc = reverse $motif1;
$motif1_rc =~ tr/ACGTNacgtn/TGCANtgcan/;
my $motif2_rc = reverse $motif2;
$motif2_rc =~ tr/ACGTNacgtn/TGCANtgcan/;

open (IN, "$FastqFile") or die "can't open $FastqFile\n";
open (OUT, ">$outfile") or die "can't open $outfile\n";
my $ReadID = "";
my $ReadSeq = "";
my $Read_qual_flag = "";
my $ReadQuality = "";

while (<IN>) {
    chomp;
    $ReadID = $_;
	$ReadSeq = <IN>;
	$Read_qual_flag = <IN>;
	$ReadQuality = <IN>;
	if ($ReadSeq =~ m/[ACGTN]*$motif1[ACGTN]+$motif2[ACGTN]*/i) {
	   print OUT ">"."$ReadID\n"."$ReadSeq\n";
    } elsif ($ReadSeq =~ m/[ACGTN]*$motif2_rc[ACGTN]+$motif1_rc[ACGTN]*/i) {
	   my $ReadSeqRC = reverse $ReadSeq;
	   $ReadSeqRC =~ tr/ACGTNacgtn/TGCANtgcan/;
	   print OUT ">"."$ReadID\n"."$ReadSeqRC\n";
	}
}
close OUT;
close IN;

system "muscle -in $outfile -out $outfile.aln.fas";
