#!/usr/bin/perl -w

# by Wenhu Guo (wenhu_guo@acgtinc.com)
#

foreach(glob("*.trimmed.rc")){
	chomp;
	$in = $_;
	$in =~ s/\.trimmed.rc//;
    system "java -Xmx8000m -jar miRanalyzer.jar input=$_ output=out/$in dbPath=/data/bin/miRanalyzer_0.3/miRanalyzerDB species=hg19 speciesShort=hsa kingdom=animal bowtiePath=/data/bin translibs=hg19_RefSeq_genes:RefSeq_genes:20:5 translibs=rfam:Rfam:20:5 translibs=repbase16:repbase16:20:5 translibs=tRNA:tRNA:20:5 translibs=pirna:pirna:20:5 reverselibs=hg19_RefSeq_genes:RefSeq_genes_rev:20:5 reverselibs=rfam:Rfam_rev:20:5 reverselibs=repbase16:repbase16_rev:20:5 reverselibs=tRNA:tRNA_rev:20:5 reverselibs=pirna:pirna_rev:20:5 minReadLength=17 minReadLengthHomolog=17 minReadLengthGenome=17 minReadLengthTrans=17 minReadLengthReverse=17 adapterTrimmed=true";
}