#!/usr/bin/perl
use strict; #use warnings;

# Created date: 02/07/2017 Samuel Moijueh
# Summarizes the taxonomy results in classification tables in each clade
# reports the results in CSV format

# Common usage: calcTaxonomy_byReads.pl Sample Sample.classify.txt 0.01

# Provide input data
die "usage:$0 <Sample.classify.txt> <percentage_threshold>\n" unless $#ARGV == 2;
my $start_run = time();

my %taxonomy_tree; # hash of hash. key1 = clade, key2 = entry, value = number of reads
my %parent_child; # stores the "parent child" relationship where the key1 is the clade, key2 is the child and the value is the parent

my %genus_parent_child; # stores the "parent child" relationship for genus rank where the key is the child and the value is the parent

my %pecentConfidence; # stores a count of the percent confidence

my $sample = shift;
my $filename = shift;
my $percentage_threshold = shift;

open my $fh, '<', $filename or die $!;
my $count = `wc -l < $filename`; chomp($count);


$taxonomy_tree{'total_number_of_roots'} = $count;
my @classification = ('domain', 'phylum', 'class', 'order', 'family', 'genus');
my %taxonomic_hierarchy = (
    "domain"  => "phylum",
    "phylum" => "class",
    "class"  => "order",
    "order"  => "family",
    "family"  => "genus",
);

my %clade_parent_child_map;

my %parent_clade_map;

my $last_entry = "";

while (<$fh>){
	chomp;
	for my $clade (@classification) {
		$_ =~  tr/"//d;
        my ($entry, $confidence) = $_ =~ /(\w+ ?\w+ ?\w+)(?:\t\Q$clade\E\t)(\S+)/i;
		my ($entry_child, $parent_clade) = get_child($clade, $_);

		$parent_clade_map{$entry} = $clade;

		if ($clade ne $parent_clade and $entry){
			$last_entry = $entry;
		}

		if ($clade ne "genus"){

			$parent_child{$clade}{$entry_child} = $entry;
		} else {
			$genus_parent_child{$entry} = $last_entry;
		}

        $taxonomy_tree{$clade}{$entry} += 1;
		push ( @{$pecentConfidence{$clade}{$entry} }, $confidence);
	}
}

sub get_child{
	my $clade = shift;
	my $line = shift;
	my ($entry_child) = $line =~ /(\w+ ?\w+ ?\w+)(?:\t\Q$taxonomic_hierarchy{$clade}\E\t)/i;
	
	if ($entry_child eq "" and $clade ne "genus"){
		return get_child($taxonomic_hierarchy{$clade}, $line);
	} else{
		return ($entry_child, $clade);
	}
}

open(my $fh_out, ">", $sample.".taxonomy_table_OTUs.csv");
print $fh_out "Total number of OTUs,$taxonomy_tree{'total_number_of_roots'}\n";

my $percentage;

# create the classification tables for the major entry in each clade
for my $clade (@classification){

	# create the headers of each table
	if ($clade eq "domain"){
		print $fh_out "Rank: Domain,Number of OTUs,Avg % Confidence, Percentage of Domain,Percent of Prior Rank, Prior Rank\n";
	} elsif ($clade eq "phylum"){
		print $fh_out "\nRank: Phylum,Number of OTUs,Avg % Confidence, Percentage of Phylum,Percent of Prior Rank, Prior Rank\n";
	} elsif ($clade eq "class"){
		print $fh_out "\nRank: Class,Number of OTUs,Avg % Confidence, Percentage of Class,Percent of Prior Rank, Prior Rank\n";
	} elsif ($clade eq "order"){
		print $fh_out "\nRank: Order,Number of OTUs,Avg % Confidence, Percentage of Order,Percent of Prior Rank, Prior Rank\n";
	} elsif ($clade eq "family"){
		print $fh_out "\nRank: Family,Number of OTUs,Avg % Confidence, Percentage of Family,Percent of Prior Rank, Prior Rank\n";
	} elsif ($clade eq "genus"){
		print $fh_out "\nRank: Genus,Number of OTUs,Avg % Confidence, Percentage of Genus,Percent of Prior Rank, Prior Rank\n";
	}

	# sort the taxonomic entries by percentages
	my @sorted_entries = sort { $taxonomy_tree{$clade}{$b} <=>
								$taxonomy_tree{$clade}{$a}} keys %{ $taxonomy_tree{$clade} };
	my $num_of_OTUs;

	for my $entry (@sorted_entries){

		if ($clade eq "domain"){
			
			$num_of_OTUs = $taxonomy_tree{$clade}{$entry};
			$percentage = $num_of_OTUs/$taxonomy_tree{'total_number_of_roots'};
			writePercentage($percentage, $entry, $num_of_OTUs, $clade, \%parent_child, \%pecentConfidence);

		} elsif ($clade eq "phylum"){
			
			$num_of_OTUs = $taxonomy_tree{$clade}{$entry};
			$percentage = $num_of_OTUs/$taxonomy_tree{'total_number_of_roots'};
			writePercentage($percentage, $entry, $num_of_OTUs, $clade, \%parent_child, \%pecentConfidence);
			
		} elsif ($clade eq "class"){

			$num_of_OTUs = $taxonomy_tree{$clade}{$entry};
			$percentage = $num_of_OTUs/$taxonomy_tree{'total_number_of_roots'};
			writePercentage($percentage, $entry, $num_of_OTUs, $clade, \%parent_child, \%pecentConfidence);

		} elsif ($clade eq "order"){
			
			$num_of_OTUs = $taxonomy_tree{$clade}{$entry};
			$percentage = $num_of_OTUs/$taxonomy_tree{'total_number_of_roots'};
			writePercentage($percentage, $entry, $num_of_OTUs, $clade, \%parent_child, \%pecentConfidence);

		} elsif ($clade eq "family"){

			$num_of_OTUs = $taxonomy_tree{$clade}{$entry};
			$percentage = $num_of_OTUs/$taxonomy_tree{'total_number_of_roots'};
			writePercentage($percentage, $entry, $num_of_OTUs, $clade, \%parent_child, \%pecentConfidence);

		} elsif ($clade eq "genus"){

			$num_of_OTUs = $taxonomy_tree{$clade}{$entry};
			$percentage = $num_of_OTUs/$taxonomy_tree{'total_number_of_roots'};
			writePercentage($percentage, $entry, $num_of_OTUs, $clade, \%parent_child, \%pecentConfidence);
		}
	}
}

# write the results to file
sub writePercentage{

	my $percentage = shift;
	my $entry = shift;
	my $num_of_OTUs = shift;
	my $clade = shift;
	my $parent_child_ref = shift;
	my $pecentConfidence_ref = shift;
	my %parent_child = %{$parent_child_ref};
	my %percentConfidence = %{$pecentConfidence_ref};

	my %clade_parent_child_map = (
	    "phylum"  => "domain",
	    "class" => "phylum",
	    "order"  => "class",
	    "family"  => "order",
	    "genus"  => "family",
	);
	

	my $upper_number;
	my $parent;
	my $sum = 0;
	my $average_confidence;

	foreach (@{$percentConfidence{$clade}{$entry}}){
		$sum += $_;
	}

	if (@{$percentConfidence{$clade}{$entry}} != 0){
		$average_confidence = $sum / @{$percentConfidence{$clade}{$entry}} * 100;
	} else {
		$average_confidence = "NA";
	}
		
	$sum = 0;

	if ($entry eq "Bacteria" || $entry eq "Archaea") {
		$upper_number = $num_of_OTUs / $taxonomy_tree{'total_number_of_roots'};
		$parent = "root";
	} else {

		if ($clade eq "genus" and $parent_child{$clade_parent_child_map{$clade}}{$entry} eq ""){
			$parent = $genus_parent_child{$entry};
		} else {
			$parent = $parent_child{$clade_parent_child_map{$clade}}{$entry};

		}

		if ($clade eq "genus" and $parent_child{$clade_parent_child_map{$clade}}{$entry} eq ""){
			$upper_number = $num_of_OTUs / $taxonomy_tree{$parent_clade_map{$parent}}{$genus_parent_child{$entry}} * 100
		} else {
			if ($taxonomy_tree{$clade_parent_child_map{$clade}}{$parent_child{$clade_parent_child_map{$clade}}{$entry}} == 0){
				if ($taxonomy_tree{$parent_clade_map{$entry}}{$parent_child{$clade_parent_child_map{$clade}}{$entry}}){
					$upper_number = $num_of_OTUs/$taxonomy_tree{$parent_clade_map{$entry}}{$parent_child{$clade_parent_child_map{$clade}}{$entry}} * 100;
				}
			} else{
				$upper_number = $num_of_OTUs/$taxonomy_tree{$clade_parent_child_map{$clade}}{$parent_child{$clade_parent_child_map{$clade}}{$entry}} * 100;

				if ($upper_number == 0) {
					$upper_number = $num_of_OTUs/$taxonomy_tree{'total_number_of_roots'};
					$parent = "root";
				}
			}
		}	
	}

	if ($percentage > $percentage_threshold){
		if ($entry ne ""){
			print $fh_out "$entry,$num_of_OTUs,$average_confidence %," . $percentage * 100 . "%," . $upper_number . "%, $parent\n";
		}
	}
}

close $fh_out;

my $end_run = time();
my $run_time = $end_run - $start_run;
print "Job took $run_time seconds\n";