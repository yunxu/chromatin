#!/bin/bash - 
for (( i=1; i <= 54; i++ )); do
  echo $i
  Rscript ENm008.chain.BoxplotWeightDist_2.R $i
done
