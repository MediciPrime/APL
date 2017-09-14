#!/usr/bin/perl -w

# by Wenhu Guo (wenhu_guo@acgtinc.com)
#


foreach(glob("*.R1.fq.gz")){
	chomp;
	$_ =~ m/(\d+_D[MN].*)\.R1\.fq\.gz/;
	system "mkdir $1";
	system "mkdir $1/fastq";
	system "mkdir $1/analysis";
	system "mv -v $1.*fq* $1/fastq";
}