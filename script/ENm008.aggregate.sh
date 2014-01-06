#!/bin/bash - 

SEQCMD="seq"
OutDir=tmp/
Cutoffs=(150 200 300 400 500 600 800 1000 1200 1500 )
#Cutoffs=(150 )

#StartSegInd=1
#EndSegInd=238
SampleSizes=(p600_d6000)
SampleSizes=($(echo p{600,900,1200}_d{3000,6000,12000,24000}))
#SeqArr=($($SEQCMD $StartSegInd $EndSegInd))

# for Cutoff in ${Cutoffs[@]}; do
# 
#   DataDir=../result/analysis/ENm008/freq/$SampleSize/cut$Cutoff
#   
#   for ((i=$StartSegInd; i<EndSegInd; i++ )); do
#       for ((j=$(($i+1)); j<=EndSegInd; j++)); do
#           ConFreq=$(grep "^$i $j 1$" $DataDir/*.con | wc -l)
#           if (( $ConFreq > 0 )); then
#               echo $i $j $ConFreq
#           fi
#       done
#   done > $OutDir/$SampleSize".cut"$Cutoff".comb.txt"
# 
# done
for SampleSize in ${SampleSizes[@]}; do
  echo $SampleSize 
  for Cutoff in ${Cutoffs[@]}; do
  
    DataDir=../result/analysis/ENm008/freq/$SampleSize/cut$Cutoff
    
    ItemID=$SampleSize.$Cutoff
    FJOB=$ItemID.q
    JOBNAME=J$FJOB
    FileLog=$0.log
  
  echo "#! /bin/bash

#$ -N $JOBNAME
#$ -S /bin/bash
#$ -o $FileLog
#$ -j y 
#$ -M yunxu
#$ -m n
#$ -cwd
 
OutFile=$OutDir/beadstring.$SampleSize.cut$Cutoff.comb.txt
TmpSortFile=\$OutFile.sorted.tmp

find $DataDir -name \*.con -exec sed /#/d {} \; | sort -k1,1n -k2,2n > \$TmpSortFile
awk 'BEGIN{SUBSEP=\"\t\"} {count[\$1,\$2] += \$3} END{for(j in count){printf (\"%s\t%.3e\n\", j,count[j])}}' \$TmpSortFile | sort -k1,1n -k2,2n > \$OutFile

  " >  $FJOB
  qsub $FJOB
  rm -f $FJOB
  
  done
done
