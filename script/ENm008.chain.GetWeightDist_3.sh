#!/bin/bash - 
# Generate contact node weighted distance

paramlist=(1_1_1_1)
OutDirRoot=/tmp
NumNodesWindow=(1 2 3 4 5 6 7 8 9)
NumNodesWindow=`seq 1 20`
NumNodesWindow=(1  3  5  8 10 15 20)
NumNodesWindow=(1 )

MaxNum=1000
for par in ${paramlist[@]}; do
  par=par.$par
  SampleSizes=`ls $par`
#  SampleSize=(K562.p1500_d3300)
  for SampleSize in ${SampleSizes[@]}; do
    OutDir=$OutDirRoot/$par/$SampleSize
    echo $OutDir
    if [ ! -d $OutDir ]; then install -d $OutDir; fi 
    for WindowSize in ${NumNodesWindow[@]}; do
      echo $WindowSize
      ItemID=$SampleSize.$WindowSize
      JOBNAME=J$ItemID.q
      FileLog=$0.log
      FJOB=$JOBNAME
(
cat <<-QSCRIPT
#$ -N $JOBNAME
#$ -S /bin/bash
#$ -o $FileLog
#$ -j y
#$ -M yunxu
#$ -m n
#$ -cwd
export PATH=$PATH:/share/apps/R/bin
Rscript  GetWeightDist_3.R $par $SampleSize $WindowSize $MaxNum 
QSCRIPT
) > $JOBNAME 

qsub $JOBNAME
rm -f $JOBNAME
    done
  done
done
