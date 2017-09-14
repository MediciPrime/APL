#!/usr/bin/perl -w

# by Wenhu Guo (wenhu_guo@acgtinc.com)
#

use Cwd;

foreach(glob("*_DM*/analysis/*.sorted.rmdup.bam")){
	chomp;
	$_ =~ m/(\d+_DM.*)\/analysis\/(\d+_DM.*)\.sorted.rmdup.bam/;
	if ("$1" ne "$2") {
	   die "Folder name and file name don't match! Fatal error.\n";
	}
	chdir "$1/analysis";
	$prefix = $1;
	$DM_sample_name = $2;
	$prefix =~ s/\_DM//;
	$DN_sample_name = "$prefix"."_DN";
	$mean_DM = `cut -f1 SV/$DM_sample_name.stats`;
	$mean_DM =~ s/\D+//;
	chomp $mean_DM;
	$stdev_DM = `cut -f2 SV/$DM_sample_name.stats`;
	$stdev_DM =~ s/\D+//;
	chomp $stdev_DM;
	$mean_DN = `cut -f1 /data/679655_LH_Leidos/Exome/Batch7/cat_set/p1/$DN_sample_name/analysis/SV/$DN_sample_name.stats`;
	$mean_DN =~ s/\D+//;
	chomp $mean_DN;
	$stdev_DN = `cut -f2 /data/679655_LH_Leidos/Exome/Batch7/cat_set/p1/$DN_sample_name/analysis/SV/$DN_sample_name.stats`;
	$stdev_DN =~ s/\D+//;
	chomp $stdev_DN;
	print "$mean_DM\t$stdev_DM\t$mean_DN\t$stdev_DN\n";
	system "/data/bin/krishna/Leidos/lumpy-sv/bin/lumpy -mw 4 -tt 0.0 -pe bam_file:$DM_sample_name.sorted.rmdup.bam,histo_file:SV/$DM_sample_name.hist,mean:$mean_DM,min_mapping_threshold:30,stdev:$stdev_DM,read_length:150,min_non_overlap:150,discordant_z:4,back_distance:20,weight:1,id:1,min_mapping_threshold:1 -pe bam_file:../../$DN_sample_name/analysis/$DN_sample_name.sorted.rmdup.bam,histo_file:../../$DN_sample_name/analysis/SV/$DN_sample_name.hist,mean:$mean_DN,stdev:$stdev_DN,read_length:150,min_non_overlap:150,discordant_z:4,back_distance:20,weight:1,id:1,min_mapping_threshold:1 > SV/$prefix.sv";
	chdir "../..";
	print "Done: $_\n";
}