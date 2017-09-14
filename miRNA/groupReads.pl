#!/usr/bin/perl -w
use strict;

	my %args; # hold the comand line arguments

# write out the header with lab information
	&header();

# parse the comand line arguments to the %args hash --> if no arguments are given --> print out the help modus
	&parseCommandLineArgs();

# the minimum number of reads of copies a sequence tag must have to be printed out
	&checkCommandLineArgs();
	## globas
	my ($outFile,$inFile,$minReadCount,$maxLen,$phred,$scarf,$eachLane,$qualT,$qualV);
	my @outputData; # this array will hold the output data
	my %count; # the count hash

	my @files = split(/,/,$inFile);
	# the count will be saved like this: $count{'read sequence'}-> array{'counts for the different samples/lanes'}  
	my $totalcount = &getCountSum(\@files,\%count);
	my @totalcount = @{$totalcount};
	
	my $nrFiles = @files;
	my @filtered;

	while(my @a = each(%count)){

#	print $a[0]." ".$a[1];
		(my $nrLanes, my $sum) = &getStat($a[1],$nrFiles,\@totalcount);
		
		if($eachLane){
			if($sum >= $minReadCount ){
				
				if($nrFiles eq $nrLanes){
					if($nrFiles > 1){
						push @outputData,"$a[0],$sum,$nrLanes";		
					}
					else{
						push @outputData,"$a[0],$sum";								
					}				
				}
			}
		}
		else{
			if($sum >= $minReadCount){
					if($nrFiles > 1){
						push @outputData,"$a[0],$sum,$nrLanes";		
					}
					else{
						push @outputData,"$a[0],$sum";								
					}				
			}
		}
	}
	
	&writeOutput($outFile,\@outputData);
	
	


###################################################################################
#######   SUBFUNCTIONS   ##########################################################
###################################################################################


# takes the count array (counts in different samples), the number of samples, and the total number of reads in each lane
sub getStat(){
	
	my @arr = @{$_[0]};
	my $sum = 0;
	my $nrLanes=0;
	for(my $i = 0; $i < $_[1]; $i++){
		if($arr[$i]){
			$sum+=$arr[$i];
			$nrLanes++;
		}
	}
	
	
	return ($nrLanes,$sum);
}

##
sub writeOutput(){
	
	open (O,">$_[0]") or die "could not open $_[0]";
	my $last = @{$_[1]};
	for(my $i = 0; $i < $last; $i++){
		my $str = $_[1]->[$i];
		$str =~ s/,/\t/g;
		print O $str."\n";
	}
	close(O);
}

sub getAllReads(){
	
	foreach my $k (keys %{$_[1]}) {
		$_[0]->{$k}++;
	}

	foreach my $k (keys %{$_[2]}) {
		$_[0]->{$k}++;		
	}
	
}


# takes an reference to an array which holds the filenames and a reference to an hash which will hold the counts
sub getCountSum(){
	
	my @totCount;
	my @fi = @{$_[0]};
	# go over all files
	for(my $i = 0; $i <= $#fi; $i++){
		print "go for $fi[$i] \n";
		my $c;
		if($scarf){
			$c = getCountsScarf($fi[$i],$_[1],$i);	
		}
		else{
			$c = getCounts($fi[$i],$_[1],$i);		
		}
		push (@totCount,$c);
	}
	return \@totCount;
}

# takes the file name, the count hash, and and index for the file 
sub getCounts(){

		if($qualT eq "mean"){
			return &hashReadsMean($_[0],$_[1],$_[2]);
		}
		elsif($qualT eq "min"){
			return &hashReadsMin($_[0],$_[1],$_[2]);			
		}
		else{
			print "something wrong with the quality input option \n\n\n";
			&commandLineOptions();
			
		}

}

sub getCountsScarf(){

		if($qualT eq "mean"){
			return &hashReadsMeanScarf($_[0],$_[1],$_[2]);
		}
		elsif($qualT eq "min"){
			return &hashReadsMinScarf($_[0],$_[1],$_[2]);			
		}
		else{
			print "something wrong with the quality input option \n\n\n";
			&commandLineOptions();
			
		}

}

sub hashReadsMeanScarf(){
	
	my $counter = 0;
	open(I,$_[0]) or die "could not open $_[0]";
	while(my $z = <I>){
		
		my @f = split(":",$z);
			
		my $seq = $f[5];
		my $qual = $f[6];
		if(length($seq) > $maxLen){
			$seq = substr($f[5],0,$maxLen);
			$qual = substr($f[6],0,$maxLen);
		}
		
		my $q = &getMeanQuality($qual);
		if($q >= $qualV){
			if($_[1]->{$seq}){
				$_[1]->{$seq}->[$_[2]]++;
			}
			else{
				my @t;
				$t[$_[2]]++;
				$_[1]->{$seq} = \@t;
			}
			$counter++;
		}
	}
	return $counter;
}


sub hashReadsMinScarf(){
	
	my $counter = 0;
	open(I,$_[0]) or die "could not open $_[0]";
	while(my $z = <I>){
		
		my @f = split(":",$z);
			
		my $seq = $f[5];
		my $qual = $f[6];
		if(length($seq) > $maxLen){
			$seq = substr($f[5],0,$maxLen);
			$qual = substr($f[6],0,$maxLen);
		}
		my $q = &getMinQuality($qual);
		print "go quality of $q\n";
		if($q >= $qualV){
			if($_[1]->{$seq}){
				$_[1]->{$seq}->[$_[2]]++;
			}
			else{
				my @t;
				$t[$_[2]]++;
				$_[1]->{$seq} = \@t;
			}
			$counter++;
		}
	}
	return $counter;
}



###########################################################
###########################################################

sub hashReadsMean(){
	
	my $counter = 0;
	open(I,$_[0]) or die "could not open $_[0]";
	while(my $z = <I>){
		
		
		$z = <I>;
		my @f = split(/\s+/," $z");
		my $seq = substr($f[1],0,$maxLen);
		
		my $header = <I>; # the quality header
		$z = <I>;
		@f = split(/\s+/," $z");
		my $qual = substr($f[1],0,$maxLen);
		
		my $q = &getMeanQuality($qual);

		if($q >= $qualV){
			if($_[1]->{$seq}){
				$_[1]->{$seq}->[$_[2]]++;
			}
			else{
				my @t;
				$t[$_[2]]++;
				$_[1]->{$seq} = \@t;
			}
			$counter++;
		}
	}
	return $counter;
}


sub hashReadsMin(){
	
	my $counter = 0;
	open(I,$_[0]) or die "could not open $_[0]";
	while(my $z = <I>){
		
		$z = <I>;
		my @f = split(/\s+/," $z");
		my $seq = substr($f[1],0,$maxLen);
#		my $seq = $f[1];
		my $header = <I>; # the quality header

		$z = <I>;
		@f = split(/\s+/," $z");
		my $qual = substr($f[1],0,$maxLen);
#		my $qual = $f[1];

		my $q = &getMinQuality($qual);
		if($q >= $qualV){
			
			if($_[1]->{$seq}){
				$_[1]->{$seq}->[$_[2]]++;
			}
			else{
				my @t;
				$t[$_[2]]++;
				$_[1]->{$seq} = \@t;
			}
			
			$counter++;
#			$_[1]->{$seq}++;
		}
	}
	return $counter;
}


#
# This function gives back the mean Q-value of the sequence/quality string
sub getMeanQuality(){
	
	my @bases = split(//,$_[0]);
	my $sum = 0;
	for(my $i = 0; $i <= $#bases; $i++){
		my $num = ord($bases[$i]) - $phred;
		$sum += $num;
	}
	
	return $sum/($#bases+1);
	
}

###
### This function gives back the Q-value of the worst base
sub getMinQuality(){
	
	my @bases = split(//,$_[0]);
	my $worst = 1000;
	for(my $i = 0; $i <= $#bases; $i++){
#		printf ("base: $bases[$i]  --> %d\n",ord($bases[$i]));
		my $num = ord($bases[$i]) - $phred;
		if($num < $worst){
			$worst = $num;
		}
	}
	return $worst;
}

sub parseCommandLineArgs(){

	if($#ARGV == -1){
		
		&commandLineOptions();
		die "\n";
	}
	for(my $i = 0; $i <= $#ARGV; $i++){
		my @f = split(/=/,$ARGV[$i]);
		$args{"$f[0]"} = $f[1];		
	}
}
	
sub header(){

print("\n\n");
print('#######################################################################'."\n");
print('###                                                                 ###'."\n");
print('###    Computational Genomics and Bioinformatics Group              ###'."\n");
print('###    Dept. of Genetics & Inst. of Biotechnology                   ###'."\n");
print('###    University of Granada, Spain                                 ###'."\n");
print('###    web: http://bioinfo2.ugr.es                                  ###'."\n");
print('###                                                                 ###'."\n");
print('###    Functional Genomics Unit, CIC bioGUNE, CIBER-HEPAD           ###'."\n");
print('###    Technology Park of Bizkaia, 48160 Derio, Bizkaia, Spain      ###'."\n");
print('###                                                                 ###'."\n");
print('###    E-mail contact: Michael Hackenberg - mlhack@gmail.com        ###'."\n");
print('###                                                                 ###'."\n");
print('###    Last update: November 2013                                   ###'."\n");
print('###                                                                 ###'."\n");
print('#######################################################################'."\n");
print("\n\n");

}	

sub commandLineOptions(){
	
	print('#######################################################################'."\n");
	print('###                                                                 ###'."\n");
	print('###  Please note that the most recent version of miRanalyzer can    ###'."\n");
	print('###  perform adapter trimming on fastq input files. We recommend    ###'."\n");
	print('###  to use miRanalyzer with a fastq input file and not with        ###'."\n");
	print('###  read/count format generated by means of this script            ###'."\n");
	print('###                                                                 ###'."\n");
	print('###  Please see also our new tool http://bioinfo5.ugr.es/sRNAbench  ###'."\n");
	print('###                                                                 ###'."\n");
	print('###  All input options are given in this format: option=value       ###'."\n");
	print('###  The following options do exist:                                ###'."\n");
	print('###                                                                 ###'."\n");
	print('###  Mandatory input paramter:                                      ###'."\n");
	print('###                                                                 ###'."\n");
	print('###  input=file1,file2,file3,...                                    ###'."\n");
	print('###  Input Files: The script accepts fastq format                   ###'."\n");
	print('###                                                                 ###'."\n");
	print('###  output=name of the output file                                 ###'."\n");
	print('###                                                                 ###'."\n");
	print('###  Optional input parameter:                                      ###'."\n");
	print('###                                                                 ###'."\n");
	print('###  minReadCount=numeric value: the minimum read count a read must ###'."\n");
	print('###  must have to be writen out (default 1)                         ###'."\n");
	print('###                                                                 ###'."\n");
	print('###  Quality values: Two different ways to control the quality      ###'."\n");
	print('###  of the reads are implemented: 1) all bases in a read must      ###'."\n");
	print('###  be better than a given threshold (Phred like Illumina scores   ###'."\n");
	print('###  --> value=10 --> one out of ten bases is incorrect, value=20   ###'."\n");
	print('###  --> one out of hundret bases is incorrect etc                  ###'."\n");
	print('###  or 2) the mean quality score must be above a given value       ###'."\n");
	print('###  The value range might be between 0 and 30                      ###'."\n");
	print('###                                                                 ###'."\n");
	print('###  the quality is set like this:                                  ###'."\n");
	print('###  qual=min:value OR                                              ###'."\n");
	print('###  qual=mean:value                                                ###'."\n");
	print('###  (default: qual=min:0)                                          ###'."\n");
	print('###                                                                 ###'."\n");
	print('###  phred=solid: the default is to rest 64 from the ASCII quality  ###'."\n");
	print('###              character (Illumina). Just set this option if the  ###'."\n");
	print('###              if the input is from solid.                        ###'."\n");
	print('###                                                                 ###'."\n");
	print('###  scarf=ture: set this option if the input is in SCARF format    ###'."\n");
	print('###                                                                 ###'."\n");
	print('###  maxLen=value: the maximal length of the reads. If set to       ###'."\n");
	print('###                higher values --> more unique reads will be      ###'."\n");
	print('###                generated (default: 26)                          ###'."\n");
	print('###                                                                 ###'."\n");
	print('###                                                                 ###'."\n");
	print('###  eachLane=true: set this option if multiple input files are     ###'."\n");
	print('###                 used. The option forces that a unique read      ###'."\n");
	print('###                 must be present in all samples (default: false) ###'."\n");
	print('###                                                                 ###'."\n");
	print('#######################################################################'."\n");
	
	
}
	
sub checkCommandLineArgs(){

	if($args{"output"}){
		$outFile = $args{"output"};
	}
	else{
		die("please specify the output with output=name");
	}
	if($args{"input"}){
		$inFile = $args{"input"};
	}
	else{
		die("please specify the output with input=name");
	}

	if($args{"minReadCount"}){
		$minReadCount = $args{"minReadCount"};
	}
	else{
		print ("Setting minimum read count to 1\n");
		$minReadCount = 1;
	}
	if($args{'maxLen'}){
		$maxLen = $args{'maxLen'}
	}
	else{
		print("Setting maximum length to 26\n");
		$maxLen = 26;
	}

	if($args{'phred'} and ($args{'phred'} eq "solid" or $args{'phred'} eq 'sanger')){
		print("Setting Phred Score type to SOLiD/Sanger (Resting 33 from ASCII)\n");
		$phred = 33;
	}
	else{
		print("Setting Phred Score type to Illumina (Resting 64 from ASCII)\n");
		$phred=64;
	}

# if scarffile
	$scarf=$args{'scarf'};
# if eachLane is set --> a read must be present in all samples (lanes)
	my $eachLane=$args{"eachLane"}; 

# the quality type 
	my $qualDef;
	if($args{'qual'}){
		$qualDef = $args{'qual'};
	}
	else{
		print("Setting quality to qual=min:0 \n");
		$qualDef="min:0";
	}
	my @t = split(/:/,$qualDef);
	$qualT = $t[0]; # the quality type
	$qualV = $t[1]; # the quality value
	
}