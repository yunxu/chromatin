#!/bin/bash - 
persistence_length_List=(600 900 1200)
nucleus_sphere_diameter_List=(6000 12000 24000)
work_path=../result/ENm008/ensemble
tmp_path=tmp
collision_length=100
packing_density=0.28
number_sample_points=640

for persistence_length in ${persistence_length_List[@]}; do
    for nucleus_sphere_diameter in ${nucleus_sphere_diameter_List[@]}; do
        keystr=p$persistence_length"_d"$nucleus_sphere_diameter
        out_path=$work_path/$keystr
        configfile=$tmp_path/ENm008.ensemble.$keystr.ini
        case $persistence_length in
          600 )
            number_nodes=238
            echo $configfile
            ;;
          900 )
            number_nodes=157
            echo $configfile
            ;;
          1200 )
            number_nodes=117
            echo $configfile
            ;;
        esac

(
cat <<-CONFFILE
  out_path = $out_path
  persistence_length = $persistence_length
  collision_length = $collision_length
  packing_density = $packing_density
  nucleus_sphere_diameter = $nucleus_sphere_diameter
  number_nodes = $number_nodes
  number_sample_points = $number_sample_points
CONFFILE
) > $configfile

#-------------------------------------------------
# begin write qsub file
#-------------------------------------------------
Prefix=($(echo {A..Z}{A..Z}))
# Ensemble size = 100000
SampleSize=1000
idx=0
while [ $idx -lt 10 ]
do
  ItemID=$keystr.${Prefix[$idx]}
  FJOB=$(basename $ItemID).q
  JOBNAME=J$FJOB.$idx
  FileLog=$0.log

(
cat <<-QSUBFILE
#! /bin/bash

#$ -N $JOBNAME
#$ -S /bin/bash
#$ -o $FileLog
#$ -j y 
#$ -M yunxu
#$ -q long.q
#$ -m n
#$ -cwd
../build/bin/ensemble -conf $configfile -prefix ${Prefix[$idx]} -samplesize $SampleSize
QSUBFILE
) > $FJOB

qsub $FJOB
rm -f $FJOB
idx=$(( $idx + 1 ))
done
#-------------------------------------------------
# end write qsub file
#-------------------------------------------------

    done
done


# Prefix=($(echo {A..Z}{A..Z}))
# #Prefix=($(echo {K..Z}{A..Z}))
# 
# # Ensemble size = 100000
# SampleSize=1000
# idx=0
# while [ $idx -lt 10 ]
# do
#   ItemID=${Prefix[$idx]}
#   FJOB=$(basename $ItemID).q
#   JOBNAME=J$FJOB.$idx
#   FileLog=$0.log
# 
# echo "
# #! /bin/bash
# 
# #$ -N $JOBNAME
# #$ -S /bin/bash
# #$ -o $FileLog
# #$ -j y 
# #$ -M yunxu
# #$ -m n
# #$ -cwd
# ../build/bin/ensemble -conf ENm008.ensemble.ini -prefix ${Prefix[$idx]} -samplesize $SampleSize
# " >  $FJOB
# qsub $FJOB
# rm -f $FJOB
# idx=$(( $idx + 1 ))
# done

