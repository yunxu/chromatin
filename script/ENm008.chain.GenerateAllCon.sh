#!/bin/bash -

PTSDIRList=(../result/ENm008/ensemble/ ../result/ENm008/chain/par.1_1_1_1/GM12878.p1500_d3300/ ../result/ENm008/chain/par.1_1_1_1/K562.p1500_d3300/)
OutDirList=(../result/analysis/ENm008/ensemble/cut.all/ ../result/analysis/ENm008/chain/cut.all/par.1_1_1_1/GM12878.p1500_d3300/ ../result/analysis/ENm008/chain/cut.all/par.1_1_1_1/K562.p1500_d3300/)

#PTSDIRList=(../result/ENm008/ensemble/p1500_d3300/) 
#OutDirList=(../result/analysis/ENm008/ensemble/cut.all/)

for (( i=0; i<${#PTSDIRList[@]}; i++)); do
	PTSDIR=${PTSDIRList[$i]}
	OutDir=${OutDirList[$i]}
  echo $OutDir
  if [ ! -d $OutDir ]; then  install -d $OutDir; fi

  ItemID=$i
  JOBNAME=J$ItemID.q
  FileLog=$0.log
  FJOB=$JOBNAME	
(	cat <<-QSCRIPT
#!/bin/bash -
#$ -N $JOBNAME
#$ -S /bin/bash
#$ -o $FileLog
#$ -j y
#$ -M yunxu
#$ -m n
#$ -cwd	
Rscript ENm008.chain.GenerateAllCon.R $PTSDIR $OutDir
QSCRIPT
) > $JOBNAME
  qsub $JOBNAME
  rm -f $JOBNAME
done
