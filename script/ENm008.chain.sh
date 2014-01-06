#!/usr/bash -
IniFiles=(ENm008.chain.GM.ini  ENm008.chain.K.ini)
Prefixs=(AA BB)


for (( Count=0; Count<${#IniFiles[@]}; Count++ )); do       
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
../build/bin/chromatin.sis.coarse -conf ${IniFiles[$Count]} -prefix ${Prefixs[$Count]}
" > $JOBNAME
qsub $JOBNAME
rm -f $JOBNAME
done

