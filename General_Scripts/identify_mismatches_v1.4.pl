#!/usr/bin/perl -w

#
# by Wenhu Guo (wenhu_guo@acgtinc.com)
# Created date: 10/20/2014
# Update date: 10/21/2014, convert a windows reference sequence file to a Unix sequence file first
#              5/1/2015, mismatches were reported based on reference position
#              5/5/2015, take one fasta sequences file instead of multiple Genbank files as input
#              5/5/2015, use grep instead of bioperl to extract query and reference sequences 
#
#

# An integrated pipeline for Nanosphere/Cepheid/SouthernResearch projects
# 1) reformat query and reference sequences;
# 2) blast query file against reference sequences;
# 3) report mismatches 

# common usage: identity_mismatches.pl query.fasta project_name reference.fasta

use strict;

die "usage $0 <query fasta file>\t<project name>\t<reference fasta file>\n" unless $#ARGV == 2;

my $query_file = shift;
my $project_name = shift;
my $ref_file = shift;

# reformat query and reference file
my $cmd_tmp = q(awk '{ sub("\r$", ""); print }');
system "$cmd_tmp $query_file > tmp";
system "fasta_formatter -i tmp -o $query_file";
system "$cmd_tmp $ref_file > tmp";
system "fasta_formatter -i tmp -o $ref_file";
system "rm -f tmp";
# make blast database
system "makeblastdb -in $ref_file -dbtype nucl -out reference";
# perform blast searches
 system "blastn -task blastn -word_size 7 -query $query_file -db reference -evalue 1e-8 -num_threads 12 > $project_name.blast";
# parse blast results and report mismatches
system "blast_to_ref_parse_WG $project_name.blast > $project_name.out";
`grep -v "Match" $project_name.out > $project_name.full_length.out`;
`grep "Match" $project_name.out > $project_name.unmatch.out`;
# produce records that need to be checked (due to mismatches at either ends of query sequences)
$cmd_tmp = q(awk 'FNR==NR{a[$1];next};!($1 in a)');
system "$cmd_tmp $project_name.full_length.out $project_name.unmatch.out > $project_name.check";

my @query_id = ();
my @query_len = ();
my @subject_id= ();
my @evalue = ();

# if the check file is not empty
if (-s "$project_name.check") {
	open (CHECKFILE, "$project_name.check") or die "Can't open $project_name.check: $!\n";
	# read in records that need to be checked and put them into arrays
	while (<CHECKFILE>) {
	  chomp;
	  my @tmp = split(/\t/,$_);
	  push @query_id, $tmp[0];
	  push @query_len, $tmp[1];
	  push @subject_id, $tmp[2];
	  push @evalue, $tmp[4];
	}
	close CHECKFILE;
    #extract query and reference sequences out if they have a match in check file
	for (my $i = 0; $i <= $#query_id; $i++) {
		`grep "$query_id[$i]" $query_file -A 1 > "$project_name.tmp_fas_$i"`;
		`grep "$subject_id[$i]" $ref_file -A 1 >> "$project_name.tmp_fas_$i"`;
	}
	my $processed_file = "$project_name.check.processed.out";
	open (PROCESSED, "> $processed_file") or die "Can't open $processed_file: $!\n";
	for (my $i = 0; $i <= $#subject_id; $i++) {
		my $tmp_file = "$project_name.tmp_fas_"."$i";
		my $aln_out = "$tmp_file".".aln.fas";
		# align query sequence and reference sequence
		system "muscle -in $tmp_file -out $aln_out -quiet";
		system "fasta_formatter -i $aln_out -o tmp";
		system "mv tmp $aln_out";
		system "rm -f $tmp_file";
		open (ALN, "$aln_out") or die "Can't open $aln_out: $!\n";
		my $identity;
		my $mismatches = "";
		# read in alignment
		while (<ALN>) {
		  chomp;
		  my $query_header = $_;
		  my $query_seq = <ALN>;
		  my $subject_header = <ALN>;
		  my $subject_seq = <ALN>;
		  my $len1 = 0;
		  my $len2 = 0;
		  my $len3 = 0;
		  if ($query_seq =~ m/([\-]{0,})([ACGTN]+[\-]{0,}[ACGTN]+)([\-]{0,})/gi) {
			  $len1 = length ($1);
			  $len2 = length ($2);
			  $len3 = length ($3);
		  }
		  my @query_seq = split(//,$query_seq);
		  my @subject_seq = split (//,$subject_seq);
		  my $count = 0;
		  for (my $i=$len1; $i <= $len1+$len2-1; $i++) {
			 if ($query_seq[$i] ne $subject_seq[$i] ) { # look for mismatches
				# position is now based on the reference sequence
				my $position = $i+1; 
				$count++;
				$mismatches .= "$query_seq[$i]>$subject_seq[$i] at position $position; ";
			 }
		 }
		 $identity = sprintf("%.2f",($len2-$count)/$len2*100); # calculate sequences identity
		}
		close ALN;
		print PROCESSED "$query_id[$i]\t$query_len[$i]\t$subject_id[$i]\t$identity\t$evalue[$i]\t$mismatches\n";
	}
	close PROCESSED;
	system "cat $project_name.full_length.out $processed_file > $project_name.final.txt";
} elsif (-z "$project_name.check") { # if the check file is empty
    system "cat $project_name.full_length.out > $project_name.final.txt";
}

print "Here is the final output:$project_name.final.txt\n";
