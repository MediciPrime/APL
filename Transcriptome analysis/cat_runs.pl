foreach(glob("run2/*_R1_001.fastq.gz")){
	chomp;
    $_ =~ m/run2\/(\d+\-R).+R1.+\.fastq\.gz/;
	$old_id = $1;
	$prefix = $old_id;
	$prefix =~ s/\-R//;
	$new_id = "$prefix"."_R";
    print "$old_id\n$new_id\n";	
	print "concatenate run1/$old_id*R1*.fastq.gz run2/$old_id*R1*.fastq.gz run3/$old_id*R1*.fastq.gz run4/$old_id*R1*.fastq.gz run5/$old_id*R1*.fastq.gz run6/$old_id*R1*.fastq.gz run7/$old_id*R1*.fastq.gz run8/$old_id*R1*.fastq.gz\n";
	system "cat run1/$old_id*R1*.fastq.gz run2/$old_id*R1*.fastq.gz run3/$old_id*R1*.fastq.gz run4/$old_id*R1*.fastq.gz run5/$old_id*R1*.fastq.gz run6/$old_id*R1*.fastq.gz run7/$old_id*R1*.fastq.gz run8/$old_id*R1*.fastq.gz > $new_id.R1.fq.gz";
	print "concatenate run1/$old_id*R2*.fastq.gz run2/$old_id*R2*.fastq.gz run3/$old_id*R2*.fastq.gz run4/$old_id*R2*.fastq.gz run5/$old_id*R2*.fastq.gz run6/$old_id*R2*.fastq.gz run7/$old_id*R2*.fastq.gz run8/$old_id*R2*.fastq.gz\n";
	system "cat run1/$old_id*R2*.fastq.gz run2/$old_id*R2*.fastq.gz run3/$old_id*R2*.fastq.gz run4/$old_id*R2*.fastq.gz run5/$old_id*R2*.fastq.gz run6/$old_id*R2*.fastq.gz run7/$old_id*R2*.fastq.gz run8/$old_id*R2*.fastq.gz > $new_id.R2.fq.gz";
}