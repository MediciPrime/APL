#!/usr/bin/perl
#
# Common usage: Remove_short_low_cov_seq.pl [filename] 1000 1
# by Wei Zhu (wei_zhu@acgtinc.com) 3/31/2016
use strict;
use warnings;

die "usage:$0 <min_length>\t<min_coverage>\n or \n $0 <file_name>\t<min_length>\t<min_coverage>\n" unless $#ARGV == 1 or $#ARGV == 2;

my $file_name = "*.scaffolds.fasta";
$file_name = shift if ($#ARGV == 2);
my $minlen = shift;
my $mincov = shift;

foreach(glob($file_name))
{  	
	open F, "<", $_;
	$_ =~ m/^(.+\.scaffolds\.)fasta/;
	my $out = $1."len_".$minlen."_cov_".$mincov.".fasta";
	print "$_	=======\>	$out\n";
	open N, ">", $out;

    local $/=">";
    while(<F>) {
        chomp;
        next unless /\w/;
        #s/>$//gs;
        my @chunk = split /\n/;
        my $header = shift @chunk;
		$header =~ /length_(\d+)_cov_(\d+)\.?/;
        my $seqlen = $1;
		my $seqcov = $2;
		#print "$header	$seqlen	$seqcov\n";

		#my $seqlen = length join "", @chunk;
        print N ">$_" if($seqlen >= $minlen && $seqcov >= $mincov);
    }
    close F;
	close N;
    local $/="\n";
}