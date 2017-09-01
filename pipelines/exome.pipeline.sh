for item in $(cat set3.list); do mv -v old/$item.* old/set3; done;

# Modify cat_run_set#.pl
perl cat_runs_set3.pl

nohup /data/bin/krishna/Leidos/FastQC/fastqc -t 24 --quiet -o FastQC *fq.gz > fastqc.nohup &

cp ../cat_set2/*pl .

parallel.sh -j 14 ../parallel_galore.pl *R1.fq.gz && /data/bin/krishna/Leidos/FastQC/fastqc -t 14 --quiet -o FastQC *val*fq.gz


mv -v ../set{2..5}/*/fastq/*_fastqc .
 
grep "Total Sequences" *R1.fq_fastqc/fastqc_data.txt | sed 's|.R1.fq_fastqc/fastqc_data.txt:Total Sequences||' >grep.raw.R1.txt
grep "Total Sequences" *R2.fq_fastqc/fastqc_data.txt | sed 's|.R2.fq_fastqc/fastqc_data.txt:Total Sequences||' >grep.raw.R2.txt
grep "Total Sequences" *R1_val_1.fq_fastqc/fastqc_data.txt | sed 's|.R1_val_1.fq_fastqc/fastqc_data.txt:Total Sequences||' >grep.trimmed.R1.txt
grep "Total Sequences" *R2_val_2.fq_fastqc/fastqc_data.txt | sed 's|.R2_val_2.fq_fastqc/fastqc_data.txt:Total Sequences||' >grep.trimmed.R2.txt

perl moving.pl 

mkdir p{1..6}
cp -a {mapping_wei.pl,VarScan_wei.pl,coverage_calc_wei.pl,flagstat_wei.pl,SV_1_wei.pl,consensus_calling_wei.pl} p1

time nohup perl mapping_wei.pl > mapping_wei.nohup &

grep "core dumped" cat_set1/Logs/*/mapping.nohup
grep "core dumped" p*/mapping_wei.nohup
grep "bam_sort_core" p*/mapping_wei.nohup
grep "Exception" p*/mapping_wei.nohup

time /data/bin/krishna/Leidos/FastQC/fastqc -t 1 --quiet *_D[MN]/fastq/*val*fq.gz &
 
time nohup perl SV_1_wei.pl > SV_1_wei.nohup &
time nohup perl flagstat_wei.pl > flagstat_wei.nohup &
time nohup perl coverage_calc_wei.pl > coverage_cal_wei.nohup &
time nohup perl VarScan_wei.pl > VarScan_wei.nohup &

time nohup perl consensus_calling_wei.pl > consensus_calling_wei.nohup &

pigz -v */fastq/*fq

# Put all samples in one folder

mv p{1..7}/*_D[MN]* .

time nohup perl cns.stat.wei.pl > cns.stat.wei.nohup &
time nohup perl variant_annotation.wei.pl > variant_annotation.wei.nohup &

# Put all paired samples in folders
# Change directories in perl scripts: Two in SV_2_wei.pl and one in CNV_wei.pl
time nohup perl SV_2_wei.pl > SV_2_wei.nohup &
time nohup perl CNV_wei.pl > CNV_wei.nohup &



# Final check & QC
mkdir Shipped Logs Logs/cns Logs/FastQC Logs/Flagstat Logs/vcf_backup

mv cat_set/p*/*D[MN] Shipped/

mv cat_set*/FastQC/*fastqc Logs/FastQC/

cd Logs/FastQC
grep "Total Sequences" *R1.fq_fastqc/fastqc_data.txt | sed 's|.R1.fq_fastqc/fastqc_data.txt:Total Sequences||' >grep.raw.R1.txt
grep "Total Sequences" *R2.fq_fastqc/fastqc_data.txt | sed 's|.R2.fq_fastqc/fastqc_data.txt:Total Sequences||' >grep.raw.R2.txt
grep "Total Sequences" *R1_val_1.fq_fastqc/fastqc_data.txt | sed 's|.R1_val_1.fq_fastqc/fastqc_data.txt:Total Sequences||' >grep.trimmed.R1.txt
grep "Total Sequences" *R2_val_2.fq_fastqc/fastqc_data.txt | sed 's|.R2_val_2.fq_fastqc/fastqc_data.txt:Total Sequences||' >grep.trimmed.R2.txt

cd ../../Shipped
grep "mapped (" */analysis/*flagstat.out > mapping_rates
grep "Unknown" */analysis/*metrics | awk '{print $1, $9}' > PCR_dup
# Remember to -24
wc -l */analysis/*vcf > variants
# Remember to /212158
wc -l */analysis/*coverage.bed > percentlocimapped

for item in $(ls */analysis/*.Illumina.exon.coverage.bed); do awk '{sum+=$NF} END {print FILENAME, sum/NR}' $item; done;

grep "%C" */analysis/*cns.stats > percent_C 
grep "%G" */analysis/*cns.stats > percent_G

cp -a */analysis/*vcf ../Logs/vcf_backup/
while read -r b a; do sed -i -e "s/Sample1/$b/" $a/analysis/$a.exome.vcf; done < name_list;
grep "Sample1" */analysis/*vcf
wc -l */analysis/*vcf > variants_after
diff variants_after variants

mv -v */analysis/*D[MN].cns.fasta* ../Logs/cns/
mv -v */analysis/*cns.stats ../Logs/cns
mv -v */analysis/*flagstat.out ../Logs/Flagstat/
rm -r -v *DM/analysis/CNV/plot

pigz */analysis/*exo*
pigz /data/679655_LH_Leidos/Exome/Batch3/Logs/cns/*fasta &

mkdir exome_batch5-{1..4}

cd exome_batch3-1
md5sum */fastq/*gz > exome_batch3-1.cksum &

rsync -av --progress exome_batch3-1 /mnt/1_1_data/

