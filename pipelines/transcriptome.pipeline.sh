mkdir p{1..8}

cp -a mRNA_analysis.pl p1

cd p1

time nohup perl mRNA_analysis.pl > p1.log &

real    13698m39.888s
user    41914m45.113s
sys     13635m50.716s

perl moving.pl

### Check 
mkdir Logs/{327..344}_R
mv *_R/*err Logs/
mv *_R/*txt Logs/
mv *_R/analysis/*.mapsplice_out/*log Logs/

grep "error" p*/*log
grep "Error" p*/*log
grep "Exception" p*/*log

fastqc -t 72 *_R/fastq/*fq.gz --quiet &

mv *_R/fastq/*_fastqc FastQC/

pigz *_R/analysis/*_R.cufflinks_out/transcripts.gtf
pigz *_R/analysis/*_R.transdecoder_out/*
pigz *_R/analysis/*bed
pigz *_R/analysis/*fa

cd Shipped/
md5sum */fastq/*gz > Transcriptome.Batch8.fastq.cksum &

wc -l */analysis/*vcf > variants


grep "Sample1" */analysis/*vcf

grep "Number of transcripts" */analysis/*transcripts.stat > transcripts

grep "unmapped_reads" */analysis/*mapsplice_out/stats.txt > mapping_rates

pigz */analysis/*vcf