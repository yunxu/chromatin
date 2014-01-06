#!/bin/bash - 

SEQCMD="seq"
OutDir=tmp/
Cutoffs=(1500)
CELLs=(GM12878 K562)
alpha=0.001
#StartSegInd=1
#EndSegInd=238
SampleSizes=(p600_d6000)
paramList=(0.1_0.1_1_1 1_0.1_1_1 0.1_10_0.1_1 )
paramList=(
    0.1_10_0.1_10
    10_0.1_0.1_10
    1_10_10_10
    1_10_1_1
  )

for param in ${paramList[@]}; do
    ParDir=par.$param
    for SampleSize in ${SampleSizes[@]}; do
      echo $SampleSize 
      for Cutoff in ${Cutoffs[@]}; do
        for CELL in ${CELLs[@]}; do
            celldir=$CELL.$SampleSizes.cut$Cutoff.a$alpha
        DataDir=../result/analysis/ENm008/chain/$ParDir/$celldir/cut$Cutoff
        
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
if [ ! -d $OutDir/$ParDir ]; then mkdir $OutDir/$ParDir; fi 
OutFile=$OutDir/$ParDir/chain.$CELL.$SampleSize.cut$Cutoff.a$alpha.comb.txt
TmpSortFile=\$OutFile.sorted.tmp

find $DataDir -name \*.con -exec sed /#/d {} \; | sort -k1,1n -k2,2n > \$TmpSortFile
awk 'BEGIN{SUBSEP=\"\t\"} {count[\$1,\$2] += \$3} END{for(j in count){printf (\"%s\t%.3e\n\", j,count[j])}}' \$TmpSortFile | sort -k1,1n -k2,2n > \$OutFile
    
      " >  $FJOB
      qsub $FJOB
      rm -f $FJOB
      
        done
      done
    done
done
