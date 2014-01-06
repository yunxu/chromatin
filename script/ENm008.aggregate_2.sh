#!/bin/bash - 

DataDir=$1
OutDir=$2
TmpSortFile=$3
TmpUnsortile=$4
OutFile=$5

#-------------------------------
#LogWeight_C=$(grep "#" $DataDir/*.con | awk '{x +=$3}END{print x/NR}' )
LogWeight_C=$(find $DataDir -name \*.con -print0 -exec grep "#" {} \; | awk '{x +=$3}END{print x/NR}' )
# echo $LogWeight_C
ExpLogWeight_C=$(echo "scale=50;e($LogWeight_C) " | bc -l)

#LogWeight_All=($(grep "#" $DataDir/*.con | cut -d" " -f 3 ))
LogWeight_All=($(find $DataDir -name \*.con -print0 -exec grep "#" {} \; | cut -d" " -f 3 ))
# LogWeight_All=($(grep "#" AA_000*.con | cut -d" " -f 3 ))

DeNorm=0
for LogWeight in ${LogWeight_All[@]}; do
    x=$(echo "e($LogWeight - ($LogWeight_C))" | bc -l)
    DeNorm=$(echo "$DeNorm + $x" | bc -l)
    # echo $DeNorm
done
DeNorm=$(echo "$DeNorm + $ExpLogWeight_C" | bc )
# echo $DeNorm

find $DataDir -name \*.con -exec grep -v "#" {} \; | sort -k1,1n -k2,2n > $TmpSortFile

perl -e '
    use strict;
    my %ContactHash;
    while(<>){
        chomp;
        my ($i, $j, $LogWeight) = split;
        # print $i." ".$j." ".$LogWeight."\n";
        my $key = $i." ".$j;
        push @{$ContactHash{$key}}, $LogWeight; 
    }
		
	my %LogWeightHash;
    foreach my $key (sort keys %ContactHash){
				$LogWeightHash{$key} = 0;
        foreach my $i (0..$#{ $ContactHash{$key}} ){
					my $x = `echo " e($ContactHash{$key}[$i] - ('$LogWeight_C'))" | bc -l`;
					$LogWeightHash{$key} += $x;
        }
				# print $LogWeightHash{$key}."\n";
				$LogWeightHash{$key} += '$ExpLogWeight_C';
				$LogWeightHash{$key} /= '$DeNorm';
				# print $LogWeightHash{$key}."\n";
    }
		
	open FILE, ">", "'$TmpUnsortile'";
	foreach my $key (sort keys %LogWeightHash){
		printf FILE "%s\t%.3e\n", $key, $LogWeightHash{$key};
	}
	close FILE;
' $TmpSortFile 

sort -k1,1n -k2,2n $TmpUnsortile > $OutFile

