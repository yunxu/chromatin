#!/bin/bash - 
scale=0.01
C4H_radii_scale=$(echo 150 "*" $scale | bc)
Radii_table_fn=radii_table


(cat <<-RADII_TABLE
"C4H" $C4H_radii_scale
RADII_TABLE
) > $Radii_table_fn

OS_NAME=$(uname)
if [ $OS_NAME == "Linux" ] ; then 
  BIN="/dump/yunxu/working/projects.new/chromatin/other/alpha/build/bin"
else
  BIN="/Users/yunxu/workspace/projects/chromatin/other/alpha/build/bin"
fi

#----------------------------------------------------------------
PTSDIR=/dump/yunxu/working/projects.new/chromatin/result/ENm008/ensemble/p1500_d3300
PTSDIR=/dump/yunxu/working/projects.new/chromatin/result/ENm008/chain/par.1_1_1_1/GM12878.p1500_d3300
#PTSDIR=/Users/yunxu/workspace/projects/chromatin/result/ENm008/chain/par.1_1_1_1/GM12878.p1500_d3300
PTSLIST=($(find $PTSDIR -name \*.pts))
MaxLogWeight=$(find $PTSDIR -name \*.pts -exec grep LogWeight {} \; | sort -k3,3n | tail -n1 | awk '{print $3}')
#MaxLogWeight=22

#----------------------------------------------------------------
probe_list=(150 270)
for PTSFN in ${PTSLIST[@]}; do
  for probe_radii in ${probe_list[@]}; do
    probe_radii_scale=$(echo $probe_radii "*" $scale | bc)
    if [[ $PTSDIR =~ "chain" ]] ; then
      BASENAME=$(basename $(dirname $PTSFN))_$(basename $PTSFN .pts).r$probe_radii
    else
      BASENAME=$(basename $PTSFN .pts).r$probe_radii
    fi

    PDBFN=$BASENAME.pdb
    python2.7 $BIN/GenPDBFromPts.py --src $PTSFN --dst $PDBFN --scale $scale --link
    $BIN/pdb2galf.pl -r $probe_radii_scale radii_table $PDBFN $BASENAME > /dev/null
    $BIN/delcx $BASENAME >/dev/null
    $BIN/mkalf $BASENAME >/dev/null
    $BIN/interface -a $BASENAME >/dev/null
    $BIN/overlap2.pl -n -l 2 $BASENAME.20.int > $BASENAME.05.a.overlap
    $BIN/postoverlap2.pl $BASENAME.05.a.overlap >/dev/null
    $BIN/volbl -s 1 $BASENAME >/dev/null
    #----------------------------------------------------------------
    LogWeight=$(grep LogWeight $PTSFN | awk '{print $3}')
    ExpDiffLogWeight=$(echo | awk -v LogWeight=$LogWeight -v MaxLogWeight=$MaxLogWeight '{printf "%.3e", exp(LogWeight-MaxLogWeight)}')
    # ExpDiffLogWeight= 3.432e-04
    #----------------------------------------------------------------
    ALFFN=$BASENAME.ovl
    echo "# ExpDiffLogWeight= $ExpDiffLogWeight" > $ALFFN
    grep -v "#" $BASENAME.05.a.overlap2 | awk '{print $1"\t"$2}' >> $ALFFN

    VOBFN=$BASENAME.vob
    echo "# ExpDiffLogWeight= $ExpDiffLogWeight" > $VOBFN
    grep -v "#" $BASENAME.1.contrib | awk '{if($8>0) print $2}' >> $VOBFN
    #----------------------------------------------------------------
    rm -f $BASENAME
    rm -f $BASENAME.1.contrib
    rm -f $BASENAME.05.a.overlap
    rm -f $BASENAME.05.a.overlap2
    rm -f $BASENAME.20.int
    rm -f $BASENAME.alf
    rm -f $BASENAME.corners
    rm -f $BASENAME.dt
    rm -f $BASENAME.info
    rm -f $BASENAME.pdb
    rm -f $BASENAME.warn

  done
#  break
done
