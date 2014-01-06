#!/usr/bin/R

rm(list=ls())

# --------------Fuction GetConFileName-----------
GetConFileName <- function(PtsFileName){
	if (length(grep("ensemble",PtsFileName))>0){
		BaseName = basename(PtsFileName)
		ConFileName = sub(".pts",".con",BaseName)
		ConFileName = paste(OutDir,"/",ConFileName,sep="")
	} else if (length(grep("chain",PtsFileName))) {
		BaseName = basename(PtsFileName)
		DirName = dirname(PtsFileName)
		ConFileName = sub(".pts",".con",paste(basename(DirName),BaseName,sep="_"))
		ConFileName = paste(OutDir,"/",ConFileName,sep="")
	}
	return(ConFileName)
}
# --------------Fuction CalPairDist-----------
CalPairDist <- function(PairVec,Pts){
	# PairVec = AllCon.Vec
	NumCon = nrow(PairVec)
	Dist = rep(0,NumCon)
	for (i in 1:NumCon){
		Dist[i] = dist(Pts[c(PairVec[i,1],PairVec[i,2]),])
	}
	return(Dist)
}

# ------------- Parameter --------------------
Cutoff = 840
# ------------- Ensemble --------------------
#PTSDIR="../result/ENm008/ensemble/"
#OutDir = "../result/analysis/ENm008/ensemble/cut.all/"
# ------------- Chain --------------------
#PTSDIR="../result/ENm008/chain/par.1_1_1_1/GM12878.p1500_d3300/"
#OutDir = "../result/analysis/ENm008/chain/par.1_1_1_1/GM12878.p1500_d3300/cut.all/"

args <- commandArgs(trailingOnly = TRUE)
PTSDIR <- args[1]
OutDir <- args[2]


PTSFILES <- list.files(path=PTSDIR,pattern=".pts",recursive=TRUE,include.dirs =TRUE)
PTSFILES <- paste(PTSDIR,PTSFILES,sep="")

#PTSFILES <- PTSFILES[1:1]
NumFiles <- length(PTSFILES)

cat("#ptsfile", length(PTSFILES),"\n")
cat("ptsfile = ", PTSFILES[1],"... \n")

# --------------Get Weight-------------------

LogWeight = rep(0,NumFiles)
for (i in 1:NumFiles){
	PtsFileName = PTSFILES[i]
	LogWeight[i] = read.table(PtsFileName,comment.char="",nrows=1)[,3]
}
MaxLogWeight = max(LogWeight)
Weight = exp(LogWeight - MaxLogWeight)

# --------------Get All Connection-----------
SampleData = read.table(PTSFILES[1])
NumNodes = nrow(SampleData)
NumPairs = (NumNodes-2)*(NumNodes-1)/2
AllCon.Vec = matrix(0,nrow=NumPairs,ncol=2)
k = 1
for (i in 3:NumNodes){
	for (j in 1:(i-2)){
		AllCon.Vec[k,] = c(j,i)
		k = k+1
	}
}

# --------------Get Distance-----------------
# --------------Inside Cutoff-----------------
for (i in 1:NumFiles){
	cat(i,"/",NumFiles,"\n")
	flush.console()
	PtsFileName = PTSFILES[i]	
	ConFileName = GetConFileName(PtsFileName)
	
	Pts = read.table(PtsFileName)
	Dist = CalPairDist(AllCon.Vec,Pts)
	InCutoffInd = which(Dist<Cutoff)

	FN = ConFileName
	# FN = "/tmp/t.txt"
	#cat(FN,"\n")
	write.table(file=FN,sprintf("# ExpDiffLogWeight= %.3e",Weight[i]),
		quote=F,row.names=F,col.names=F,append=F)
	write.table(file=FN,AllCon.Vec[InCutoffInd,],quote=F,row.names=F,col.names=F,append=T)
}

