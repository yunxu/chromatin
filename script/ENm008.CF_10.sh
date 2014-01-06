#!/usr/bash -
SampleSizeList=($(echo p{1500,600}_d{3300,4400,6600,8800}))
SampleSizeList=($(echo p1500_d3300))
for SampleSize in ${SampleSizeList[@]}; do
    ItemID=$SampleSize
    JOBNAME=J$ItemID.q
    FileLog=$0.log
    FJOB=$JOBNAME
    echo "
#$ -N $JOBNAME
#$ -S /bin/bash
#$ -o $FileLog
#$ -j y
#$ -M yunxu
#$ -m n
###$ -q long.q          
#$ -cwd
export PATH=$PATH:/share/apps/R/bin
Rscript ENm008.CF_10.R $SampleSize
    " > $JOBNAME
    qsub $JOBNAME
    rm -f $JOBNAME
done
