#!/bin/bash - 

CELLList=(GM12878 K562)
SampleSizeList=(p1500_d{3300,4400,6600,8800})

TwoLetterList=($(echo {A..Z}{A..Z}{A..Z}))
TwoLetterSubList=(${TwoLetterList[@]:0:1000})

number_sample_points=360
OutRootDir=/dump/temp
OutRootDir=../result/ENm008/chain/
#OutRootDir=/dump/yunxu/working/projects.new/chromatin/result/ENm008/chain

M_max=100

paramList=($(echo {0.1,1,10}_{0.1,1,10}_{0.1,1,10}_{0.1,1,10}))
paramList=(
    0.1_10_0.1_10
    10_0.1_0.1_10
    1_10_10_10
    1_10_1_1
      )
Count=0

for param in ${paramList[@]}; do
  rho_1=$(echo $param | awk -F"_" '{print $1}')
  rho_2=$(echo $param | awk -F"_" '{print $2}')
  rho_3=$(echo $param | awk -F"_" '{print $3}')
  tau_t=$(echo $param | awk -F"_" '{print $4}')
  param="par."$param
    for SampleSize in ${SampleSizeList[@]}; do
      persistence_length=$(echo $SampleSize | sed -e 's/p\(.*\)_d\(.*\)/\1/')
      nucleus_sphere_diameter=$(echo $SampleSize | sed -e 's/p\(.*\)_d\(.*\)/\2/')
      for CELL in ${CELLList[@]}; do
                KeyBase=$CELL.$SampleSize
                ConfFile=$KeyBase.$param.conf
                OutDir=$OutRootDir/$param/$KeyBase
                if [ ! -d $OutDir ]; then install -d $OutDir; fi
                echo "
                out_path = $OutDir
                persistence_length = $persistence_length 
                collision_length = 300
                packing_density = 0.28
                cutoff = 500
                number_sample_points = $number_sample_points
                nucleus_sphere_diameter = $nucleus_sphere_diameter
                start_end_file = data/analysis/ENm008/ENm008.p$persistence_length.equ.seg.len.txt
                pval_file = data/analysis/ENm008/$KeyBase.q5.a5.sis.seg.dst
                M_max = $M_max
                rho_1 = $rho_1
                rho_2 = $rho_2
                rho_3 = $rho_3
                tau_t = $tau_t
                " > $ConfFile
                continue
              for TwoLetter in ${TwoLetterSubList[@]}; do            
    
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
##$ -q short.q
#$ -cwd
../build/bin/chromatin.sis.coarse -conf $ConfFile -prefix $TwoLetter
                          " > $JOBNAME
                          qsub $JOBNAME
                          rm -f $JOBNAME
                
                          done
                  done
                done

done

