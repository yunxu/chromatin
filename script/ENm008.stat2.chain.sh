#!/usr/bash -

Cutoffs=(1500 )
PerlCommand="ENm008.stat2.chain.pl"
SampleSizes=(p600_d6000)
CELLs=(GM12878 K562)
alpha=0.001
paramList=(0.1_0.1_1_1 1_0.1_1_1 0.1_10_0.1_1 )
paramList=(
    0.1_10_0.1_10
    10_0.1_0.1_10
    1_10_10_10
    1_10_1_1
  )


for param in ${paramList[@]}; do
pardir=par.$param
    for SampleSize in ${SampleSizes[@]}; do
      for CELL in ${CELLs[@]}; do
        case $param in
          0.1_0.1_1_1 )
            case $CELL in
                # quantile 25%
               GM12878 )
                 ErrCut=3.079905
                 ;;
               K562 )
                 ErrCut=2.315960
                 ;;
            esac
            ;;
          0.1_10_0.1_1 )
            case $CELL in
                # quantile 25%
               GM12878 )
                 ErrCut=2.351527
                 ;;
               K562 )
                 ErrCut=1.352552
                 ;;
            esac
            ;;
          1_0.1_1_1 )
            case $CELL in
                # quantile 25%
               GM12878 )
                 ErrCut=2.960933
                 ;;
               K562 )
                 ErrCut=2.158715
                 ;;
            esac
            ;;
        esac  
       echo $ErrCut" "$pardir" "$CELL 
          for Cutoff in ${Cutoffs[@]}; do
        celldir=$CELL.$SampleSize.cut$Cutoff.a$alpha
        
        EnsemblePath="../result/ENm008/chain/$pardir/$celldir"
        TmpFile=/tmp/$pardir.$celldir.tmp
        find $EnsemblePath -name 0001.pts -exec echo -n {} \; -exec tail -n2 {} \; | paste -s -d" \n"  > $TmpFile
        MaxLogWeight=$(awk '{print $3}' $TmpFile | sort -k1,1n | tail -n1)
#        MaxLogWeight=$(find $EnsemblePath -name 0001.pts -print0 -exec grep "# LogWeight" {} \;| cut -d" " -f 3 | sort -n | tail -n1)
      OutPath=../result/analysis/ENm008/chain/$SampleSize
      if [ ! -d $OutPath ]; then install -d $OutPath;  fi
#      lsfiles=($(ls $EnsemblePath))
        lsfiles=($(awk '{print $1}' $TmpFile| xargs -n1 dirname | xargs -n1 basename ))
        for Prefix in ${lsfiles[@]};
        do
        	echo $Prefix
          # if [ ! -d $OutPath/$Prefix ]; then mkdir $OutPath/$Prefix; fi
          ItemID=$Prefix
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
#$ -q long.q
perl $PerlCommand $SampleSize $Prefix $Cutoff $MaxLogWeight $pardir $celldir
        " > $JOBNAME
          qsub $JOBNAME
          rm -f $JOBNAME
    #    break
          done
        done
      done
    done
done
  
