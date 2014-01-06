#!/bin/bash - 

SEQCMD="seq"
OutDir=tmp/
Cutoffs=(150 200 300 400 500 600)
#Cutoffs=(150 )

#StartSegInd=1
#EndSegInd=238
SampleSizes=(p600_d6000)
SampleSizes=($(echo p{600,900,1200}_d{6000,12000,24000}))
SampleSizes=($(echo p{600,900}_d{6000,12000,24000}))
#SeqArr=($($SEQCMD $StartSegInd $EndSegInd))

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

DataDir=$DataDir
OutDir=$OutDir 
OutFile=$OutDir/beadstring.$SampleSize.cut$Cutoff.comb.txt
TmpSortFile=\$OutFile.sorted.tmp
TmpUnsortFile=\$OutFile.unsorted.tmp

sh ENm008.aggregate_2.sh $DataDir $OutDir \$TmpSortFile \$TmpUnsortFile \$OutFile
  " >  $FJOB
  qsub $FJOB
  rm -f $FJOB
  done
done
