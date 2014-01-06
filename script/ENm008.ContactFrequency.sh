#!/usr/bash -

NumIDs=($(seq -f%03g 1 100))

for i in ${NumIDs[@]}; do
    echo $i
    FileLog=ENm008.ContactFrequency.sh.log
    JOBNAME=J$i
    QName=$JOBNAME.q
    echo "
#$ -N $JOBNAME
#$ -S /bin/bash
#$ -o $FileLog
#$ -j y
#$ -M yunxu
#$ -m n
#$ -cwd

Rscript ENm008.ContactFrequency.R $i
" > $QName
qsub $QName
rm -f  $QName
done
