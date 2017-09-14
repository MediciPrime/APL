#!/usr/bin/perl

=head1 AUTHOR

 Dimitar Kenanov, <dimitark@aitbiotech.com>
 Modified by Wenhu Guo, <wenhu_guo@acgtinc.com> on 06/28/2014

=cut

use strict;
use warnings;
use Getopt::Long;
use List::Util qw(max);
use Pod::Usage;

my ($opts,$help,$man,$bam_in,$trans_in,$trans_out,$pile_out,$err_file,@chroms);

$opts=GetOptions("b=s"=>\$bam_in,
				 "t=s"=>\$trans_in,
				 "h"=>\$help) or pod2usage(2);

pod2usage(1) if($help);
pod2usage("$0: No files given!\n")  if ((!$bam_in || !$trans_in));

$trans_out="transcripts.fa";
$pile_out="$bam_in.v7.mpileup";

# file to show transcript with no start or end
$err_file="$bam_in.v7.log";

my ($fhi,$fho,$fhe,$line,$region,$gid,$rid,$fpkm,$nreads,$found_seqs,$ntrans,@tmp);

open($fho,'>',$trans_out) or die "problems $trans_out: $!\n";
open($fhe,'>',$err_file) or die "problems $err_file: $!\n";

$|=1;

$found_seqs=0;
$ntrans=0;

open($fhi,'<',$trans_in) or die "problems $trans_in: $!\n";
while($line=<$fhi>){
	chomp $line;
	$ntrans++;
	@tmp=split/\t/,$line;
	$rid=$tmp[3];
	$gid=$tmp[8];
	my $seq;
	if ($tmp[9] > 0) {
	    $seq=extract_seq($tmp[0],$tmp[1],$tmp[2]);
		print $fho ">$rid,$gid,$tmp[0]:$tmp[1]-$tmp[2],$tmp[9]reads\n$seq\n";
		$found_seqs++;
	}else{
		print $fhe "NO SEQ FOUND FOR REGION: $rid,$gid,$tmp[0]:$tmp[1]-$tmp[2],$tmp[9]reads\n";
	}
	print "\n";
}
close $fhi;

print "VALID TRANSCRIPTS: $ntrans :\n";
print "FOUND SEQS: $found_seqs :\n";

close $fho;
close $fhe;

###################
###### SUBS #######
###################

sub extract_seq{
	my ($chr,$start,$end)=@_;
	my ($region,$res,$ret_seq,$n,$i);
	my (@seq,@pileup);
	$region="$chr:$start-$end";
	print "$region\n";
	$res=`samtools_ndsm mpileup -ABQ0 -d 1000000 -r $region $bam_in> $pile_out`;	
	put_in_memory($pile_out,\@pileup,$chr);
	foreach $i (@pileup){
		process_line(\@$i,\@seq,\$n);
	}
	$ret_seq=join('',@seq);
	return $ret_seq;
}

sub put_in_memory{
	my ($file,$arr,$chr)=@_;
	my ($fh,$line,@tmp);
	open($fh,'<',$file) or die "problems $file:$!\n";
	while($line=<$fh>){
		chomp($line);
		@tmp=split/\t/,$line;
		if($tmp[0] eq $chr){
			push @$arr,[@tmp[0..4]];
		}
	}
	close $fh;
}

sub process_line{
	my ($inarr,$arr,$n)=@_;
	my @tmp;
	my ($sa,$sc,$sg,$st);
	my (@a,@c,@g,@t);
	my %svals;
	@a=$$inarr[4] =~/(a)/ig;
	@c=$$inarr[4] =~/(c)/ig;
	@g=$$inarr[4] =~/(g)/ig;
	@t=$$inarr[4] =~/(t)/ig;
	$sa=scalar(@a);
	$sc=scalar(@c);
	$sg=scalar(@g);
	$st=scalar(@t);
	$svals{'a'}=$sa;
	$svals{'c'}=$sc;
	$svals{'g'}=$sg;
	$svals{'t'}=$st;
	
	my $max=0;
	my $letter='';
	my $k;
	foreach $k (keys %svals){
		if($svals{$k}>$max){
			$max=$svals{$k};
			$letter=$k;
		}
	}		
	
	if($max && $letter){
		push @$arr,uc($letter);
	}else{
		$$n++;
	}
}



