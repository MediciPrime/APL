#!/usr/bin/perl -w
use strict;

# Originally from JPM, modified and updated by WG
# Pull out a subset of paired-end reads by fraction
# Version update:
# 1.0, 2/16/2015

die "usage: $0 <R1_reads> <R2_reads> <Project_name> <Fraction-Reads-To-Keep:0-1>\n" unless $#ARGV == 3;
my $R1_reads = shift;
my $R2_reads = shift;
my $Project_name = shift;
my $FracToKeep = shift;

my $KeepReads = 0;
my $LineCountMax = 4;
my $LineCount = $LineCountMax;

open (FILEA, "< $R1_reads") or die "can't open $R1_reads\n";
open (FILEB, "< $R2_reads") or die "can't open $R2_reads\n";
open (OUTFILEA, "> $Project_name.R1.$FracToKeep.fq") or die "can't open $Project_name.R1.$FracToKeep.fq\n";
open (OUTFILEB, "> $Project_name.R2.$FracToKeep.fq") or die "can't open $Project_name.R2.$FracToKeep.fq\n";

while (!eof(FILEA) and !eof(FILEB)) {
    my $lineA = <FILEA>;
	my $lineB = <FILEB>;
	
    if ($LineCount == $LineCountMax) {
	$LineCount = 0;
	$KeepReads = 0;
	$KeepReads = 1 if rand() < $FracToKeep;
    }
	
    if ($KeepReads) {
	print OUTFILEA $lineA;
	print OUTFILEB $lineB; 
    }

    $LineCount++;
}

close FILEA;
close FILEB;
close OUTFILEA;
close OUTFILEB;
