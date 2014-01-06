#!/usr/bin/R
rm(list=ls())
#------------------------------------------------------------------
library(grid)
library(gridExtra)
library(lattice)
library(RColorBrewer)
library(reshape)
library(ggplot2)
#------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (Sys.info()["sysname"] == "Darwin"){
  CELL = "K562"
  # CELL = "GM12878"
}else{
  CELL = args[1]
}
OutDir = "tmp"
#------------------------------------------------------------------
# NumMax = 1000
Type="r270.ovl"
method="alpha"
#------------------------------------------------------------------
Type="con"
method = "cut.all"
#------------------------------------------------------------------
ConDir=paste("/Users/yunxu/workspace/projects/chromatin/result/analysis/ENm008/chain/",method,"/par.1_1_1_1/",CELL,".p1500_d3300",sep="")
# ConDir=paste("/dump/yunxu/working/projects.new/chromatin/result/analysis/ENm008/chain/",method,"/par.1_1_1_1/",CELL,".p1500_d3300",sep="")
SegLenFN="data/analysis/ENm008/ENm008.p1500.equ.seg.len.txt"
NumNodes=nrow(read.table(SegLenFN))+1
#------------------------------------------------------------------
GetWeightFile <- function(ConDir){
  ConFiles   = list.files(path=ConDir,pattern=Type,full.names=T)
  NumFiles = length(ConFiles)
  # NumFiles = 100
  Weight = rep(0,NumFiles)
  ConList = list()
  for (i in 1:NumFiles){
    Weight[i] = read.table(ConFiles[i],comment.char="",nrows=1)[,3]
    ConList[[i]] = read.table(ConFiles[i])
  }
  return(list("weight"=Weight,"conlist"=ConList))
}
#------------------------------------------------------------------
#------------------------------------------------------------------

GetPIJMat <- function(ConDir){
  WeightData <- GetWeightFile(ConDir)
  Weight <- WeightData$weight
  ConList <- WeightData$conlist
  NumFiles <- length(Weight)
  #------------------------------------------------------------------
  NewPIJMat = matrix(NA,nrow=NumNodes,ncol=NumNodes)
  for (i in 1:(NumNodes-2)){
    for (j in (i+2):NumNodes){
      cat(i, " / ",j ," / ", NumNodes,"\n"); flush.console();
      ConIJ = unlist(lapply(ConList,function(x){length(which(x[,1]==i & x[,2]==j))}))
      NewPIJMat[i,j] = sum(ConIJ * Weight) / sum(Weight)
    }
  }
  return(NewPIJMat)
}

#------------------------------------------------------------------
write.alltriplet.matrix <- function(m, file="", append=F, fileEncoding=""){
  if (file == ""){
    file <- stdout()
  }
  write.table(paste("#",nrow(m),ncol(m)),file,append=F,quote=F,row.names=F,col.names=F)
  for (i in 1:(nrow(m)-1)){
    for (j in (i+1):ncol(m)){
      if(!is.na(m[i,j]) ){
        write.table(sprintf("%d\t%d\t%.3e",i,j,m[i,j]),file,append=T,quote=F,row.names=F,col.names=F)
      }
    }
  }
}
#------------------------------------------------------------------
read.alltriplet.matrix <- function(file){
  # file = FN
  HeadInfo = read.table(file,comment.char="",nrows=1)
  NumRow = HeadInfo[1,2]
  NumCol = HeadInfo[1,3]

  m = matrix(0,nrow=NumRow,ncol=NumCol)
  DF = read.table(file)
  for (i in 1:nrow(DF)){
    m[DF[i,1],DF[i,2]] = DF[i,3]
    m[DF[i,2],DF[i,1]] = DF[i,3]
  }

  return(m)
}
#------------------------------------------------------------------

(FN <- paste(OutDir,"/",CELL,".p1500_d3300.c840.",Type,".pij.triplet.allall.txt",sep=""))
###################################################################
# NewPIJMat <- GetPIJMat(ConDir)
# write.alltriplet.matrix(NewPIJMat,file=FN)
###################################################################
NewPIJMat <- read.alltriplet.matrix (FN)
###################################################################
translate.pij2dist <- function(pijmat){
  MaxDist = 3300
  MinDist = 300
  SegLen = read.table(SegLenFN)

  DistMat = matrix(0,nrow=nrow(pijmat),ncol=ncol(pijmat))
  for (i in 1:(nrow(pijmat)-1)){
    for (j in (i+1):(ncol(pijmat))){
      if ( j - i == 1){
        DistMat[i,j] = DistMat[j,i] = SegLen[i,1]
      }else{
        DistMat[i,j] = DistMat[j,i] = MaxDist + (MinDist - MaxDist) *
        pijmat[i,j]
      }
    }
  }
  return(DistMat)
}
#------------------------------------------------------------------
Node.Dist <- as.dist(translate.pij2dist(NewPIJMat))
mds1 = cmdscale(Node.Dist, k=3)


plot(mds1[,1], mds1[,2], type = "n", xlab = "", ylab = "", axes = FALSE,
       main = "cmdscale (stats)")
text(mds1[,1], mds1[,2], 1:nrow(mds1), cex=0.9)

