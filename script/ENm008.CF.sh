#!/usr/bash -
CELLList=(K562 GM12878)
#SampleSizeList=(640_r1100 640_r6000 640_r60000)
SampleSizeList=(640_r6000)
#SampleSizeList=(640_r6000 )
alphaList=(0.5)

Count=0
for CELL in ${CELLList[@]}; do
  for SampleSize in ${SampleSizeList[@]}; do
      for alpha in ${alphaList[@]}; do
       Count=$(($Count+1))
  ItemID=$Count
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
#$ -cwd
Rscript ENm008.CF.R $CELL $SampleSize $alpha
" > $JOBNAME
qsub $JOBNAME
rm -f $JOBNAME
      done
  done
done

