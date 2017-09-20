#!/usr/bin/perl -w
## Created date: 01/29/2016 by Wei Zhu
## Modified: 02/13/2017 by Wei Zhu: fixed a format issue (Sample OTUs have "-" ahead of "Root") that would confuse the Krona 

use strict;

die "usage:$0 <project.classify.txt>\t<project.uc>\n" unless $#ARGV == 2;

my $project_name = shift;
my $file1 = shift;
my $file2 = shift;

my $count = 0;

open FILE1, "<", $file1;

open FILE2, "<", $file2;

open FILE3, ">", $project_name.".classify_uc.txt";

my %hash = ();

while (my $line1 = <FILE1>){
	my ($key, $value) = (split(/\s+/, $line1, 2))[0, 1]; 
	#print "$key	$value\n";
	$hash{ $key } = $value; 
} 



while (my $line2 = <FILE2>){
	my ($read, $otu) = (split(/\s+/, $line2))[8, 9];

	#print "$read	$otu\n";
	$count++;
	
	if ($hash{ $otu}=~ /-\tRoot/)
		{print FILE3 "$read	$hash{ $otu}";}
	else 
		{print FILE3 "$read		$hash{ $otu}";}

}

print "$count\n";


close FILE1;
close FILE2;
close FILE3;