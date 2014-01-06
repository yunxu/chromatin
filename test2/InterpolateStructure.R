#!/usr/bin/R

SOURCE = "GM.pts"
TARGET = "K.pts"
OutDir = "pts"
s.pts = read.table(SOURCE)[,1:3]
t.pts = read.table(TARGET)[,1:3]

Nstep = 100
# do.breaks(c(s.pts[5,1],t.pts[5,1]),Nstep)

Arr = array(0,c(nrow(s.pts),3,Nstep))

for (iRow in 1:nrow(s.pts)){
	for (j in 1:3){
		Arr[iRow,j,] = do.breaks(c(s.pts[iRow,j],t.pts[iRow,j]),Nstep-1)
	}
}

for (i in 1:Nstep){
	FN = sprintf("%s/mov_%03.0f.pts", OutDir,i)
	write.table(file=FN,Arr[,,i],quote=FALSE,row.names=F,col.names=F)
}



Arr = array(0,c(nrow(s.pts),3,Nstep))

for (iRow in 1:nrow(s.pts)){
	for (j in 1:3){
		Arr[iRow,j,] = do.breaks(c(t.pts[iRow,j],s.pts[iRow,j]),Nstep-1)
	}
}
for (i in 1:Nstep){
	FN = sprintf("%s/mov_%03.0f.pts", OutDir,i+Nstep)
	write.table(file=FN,Arr[,,i],quote=FALSE,row.names=F,col.names=F)
}


#!/bin/bash
# for file in pts/*; do pdbfile=pts/$(basename $file .pts).pdb; python GenPDBFromPts.py $file $pdbfile; done
# find pts -name \*.pdb -exec echo "load " {} ", mov" \; > /tmp/list

bg_color white
set sphere_scale,10
spectrum
select hs40, resi 14
select hs10, resi 17
select a2, resi 20
show sphere,  a2 or hs40 or hs10
color red, a2
color pink, hs10
color blue, hs40

set ray_trace_frames, 1 
show sticks
set stick_radius, 3