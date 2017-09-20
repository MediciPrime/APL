#!/usr/bin/perl -w

#
# Modified from "metabarcoding.commands_WG.txt" and
# "resequencing_pipeline_v1.5_ngs2.pl"
# 
# Created date: 01/06/2016 Wei Zhu
# Updated: 10/13/2016 Wei Zhu v1.2
# Updated: 03/02/2017 Samuel Moijueh v1.3
# An integrated pipeline for metabarcoding analysis
# 1) trim raw reads if not already trimmed (trim_galore and/or sickle) 
# 2) merge FastQ reads (pandaseq)
# 3) quality-filter and convert FastQ  to Fasta (usearch)  
# 4) dereplicate Fasta reads (vsearch)
# 5) cluster otus (usearch)
# 6) filter chimera (usearch) 
# 7) map otus to original files (vsearch)
# 8) classify otus (classify.jar)
# 9) generate otu table (uc2otutab.py)
# 10)provide HTML file classify, hierarchy, and otus (ktImportRDP)
# 11) provide tabular taxonomic classification based on OTUs or reads(calcTaxonomy_byOTUs_v1.0.pl and calcTaxonomy_byReads_v1.0.pl)

# Common usage: metabarcoding.pl project_name sample.R1.fq sample.R2.fq 1 1

use strict;

# Provide project name and input data
die "usage:$0 <project_name>\t<read1_forward>\t<read2_reverse>\t<trim_galore_trimming: 0 or 1>\t<sickle_trimming: 0 or 1>\n" unless $#ARGV == 4;

my $project_name = shift;
my $ref_db = "/data/reference/fungalits_UNITE_trainingdata_07042014/UNITE.RDP_04.07.14.rmdup.fasta";
my $read1 = shift;
my $read2 = shift;

#my $forward_primer = shift;
#my $reverse_primer = shift;
#my $threads = shift;
my $trim_galore_trimming = shift;
my $sickle_trimming = shift;

my $pandaseq_merge = 1;
my $read1_trimmed;
my $read2_trimmed;
my $single_reads;

# if reads were already trimmed
unless ($trim_galore_trimming || $sickle_trimming) {
       $read1_trimmed = $read1;
	   $read2_trimmed = $read2;
}

# trim reads with trim_galore
if ($trim_galore_trimming) {
   print "***Trim_galore_trimming running...\n";	
   system "trim_galore -q 20 --length 100 --no_report_file --paired $read1 $read2 2>$project_name.trim_galore.err";
   if ($read1 =~ /\.fastq$/ && $read2 =~ /\.fastq$/){
       $read1 =~ s/\.fastq$/_val_1.fq/;
	   $read2 =~ s/\.fastq$/_val_2.fq/;
   } elsif ($read1 =~ /\.fq$/ && $read2 =~ /\.fq$/){
       $read1 =~ s/\.fq$/_val_1.fq/;
	   $read2 =~ s/\.fq$/_val_2.fq/;
   } else{
       $read1 =~ s/\.fastq$/_val_1.fq/;
	   $read2 =~ s/\.fastq$/_val_2.fq/;
   }
   $read1_trimmed = $read1;
   $read2_trimmed = $read2;
}

# trim reads with sickle (necessary if raw reads have low quality scores)
if ($sickle_trimming) {
   if ($trim_galore_trimming) {
	print "***sickle_trimming running...\n";
      	my $prefix1 = $read1_trimmed;
	  my $prefix2 = $read2_trimmed;
	  $prefix1 =~ s/_val_1\.fq//;
	  $prefix2 =~ s/_val_2\.fq//;
	  $read1_trimmed = $prefix1."_val_1.sickle.fq";
	  $read2_trimmed = $prefix2."_val_2.sickle.fq";
	  $single_reads = $prefix1.".sickle.single.fq";
      system "sickle pe -q 20 -l 100 -t sanger -f $read1 -r $read2 -o $read1_trimmed -p $read2_trimmed -s $single_reads 1>$project_name.sickle.log";
	}
}

# merge reads with pandaseq
# Forward	ITS-1	TCCGTAGGTGAACCTGCGG
# Forward	ITS-5	GGAAGTAAAAGTCGTAACAAGG
# Reverse	ITS-4	TCCTCCGCTTATTGATATGC
# Forward	NL-1	GCATATCAATAAGCGGAGGAAAAG
# reverse	NL-4	GGTCCGTGTTTCAAGACGG
# Forward	ITS1.BITS	ACCTGCGGARGGATCA
# reverse	ITS1.B58S3	GAGATCCRTTGYTRAAAGTT

if ($pandaseq_merge) {
     print "\n***pandaseq_merge running...\n\n";	
     system "pandaseq -F -N -T 10 -l 50 -o 10 -f $read1_trimmed -r $read2_trimmed -w $project_name.merged.fastq -g $project_name.pandaseq.log";
}

# quality filtering and fastq-fasta
# pick a max error rate (maxee) between 0.2 and 0.5
print "\n***usearch -fastq_filter running...\n\n";
system "usearch8.1.1861_i86linux32 -fastq_filter $project_name.merged.fastq -fastaout $project_name.merged.fasta -fastq_maxee_rate 0.2";

# add barcode label to each fasta file
# THIS STEP IS ESSENTIAL
system "sed -i 's/>.*/&;barcodelabel=$project_name;/' $project_name.merged.fasta";
    
#Concatenate all merged fasta files from same run
#system "cat *.merged.fasta > $project_name.merged.fasta" 

#if pooled ITS and 16S, split using first 

# dereplicate
print "\n***vsearch --derep_fulllength running...\n\n";
system "vsearch --derep_fulllength $project_name.merged.fasta --output $project_name.derep.fasta -sizeout -minseqlength 64"; 

# cluster otus
print "\n***usearch -cluster_otus running...\n\n";
system "usearch8.1.1861_i86linux32 -cluster_otus $project_name.derep.fasta -otus $project_name.otus.fasta -minsize 2 -relabel OTU_ -sizeout -uparseout $project_name.UPARSE.cluster.txt";

# chimera filter- use latest RDP training set from http://sourceforge.net/projects/rdp-classifier/files/RDP_Classifier_TrainingData/
# use RDPClassifier_16S_trainset for 16S
# use fungalits_UNITE_trainingdata for ITS
print "\n***vsearch --uchime_ref running...\n\n";
system "vsearch --uchime_ref $project_name.otus.fasta --db $ref_db --strand plus --minh 1.0 --nonchimeras $project_name.otus.good.fasta --uchimeout $project_name.otus.uchime --uchimealns $project_name.otus.aln";

# map otus to original files
print "\n***vsearch --usearch_global running...\n\n";
system "vsearch --usearch_global $project_name.merged.fasta -db $project_name.otus.good.fasta -id 0.97 --alnout $project_name.alnout.txt --uc $project_name.uc";

# generate otu table
system "python2.7 /data/bin/uc2otutab.py $project_name.uc > $project_name.otu.tbl.txt";

# classify otus (RDP Classifier Version 2.12)
print "\n***classifier.jar classify running...\n\n";
system "java -Xmx1g -jar /data/bin/rdp_classifier_2.12/dist/classifier.jar classify -c 0.80 -g fungalits_warcup -h $project_name.hierarchy.txt -o $project_name.classify.txt $project_name.otus.good.fasta";

# provide HTML file classify, hierarchy, and otus
print "\n***ktImportRDP running...\n\n";
system "ktImportRDP -o $project_name.classify.html $project_name.classify.txt";

# add read counts to classiy.txt
system "combine_uc_classify_v1.1.pl $project_name $project_name.classify.txt $project_name.uc";

# provide HTML file classify, hierarchy, and reads
system "ktImportRDP -o $project_name.classify_uc.html $project_name.classify_uc.txt";

# provide tabular taxonomic classification based on OTUs or reads
system "calcTaxonomy_byOTUs_v1.1.pl $project_name $project_name.classify.txt 0";
system "calcTaxonomy_byReads_v1.1.pl $project_name $project_name.classify_uc.txt 0";
