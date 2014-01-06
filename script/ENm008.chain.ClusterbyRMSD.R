#!/bin/R
rm(list=ls())
#------------------------------------------------------------------
library(fpc)
args <- commandArgs(trailingOnly = TRUE)


if (Sys.info()["sysname"] == "Darwin"){
	CELL = "GM12878"
	CELL = "K562"
	DataDir = paste("/Users/yunxu/workspace/projects/chromatin/result/analysis/ENm008/chain/weight.dist/par.1_1_1_1/",CELL,".p1500_d3300",sep="")
	ContactFile = paste("/Users/yunxu/workspace/projects/chromatin/script/data/analysis/ENm008/",CELL,".p1500_d3300.c840.a5.sis.node.dst",sep="")
} else {
	CELL = args[1]
	DataDir = 
	paste("/dump/yunxu/working/projects.new/chromatin/result/analysis/ENm008/chain/weight.dist/par.1_1_1_1/",CELL,".p1500_d3300",sep="")
	ContactFile = paste("/dump/yunxu/working/projects.new/chromatin/script/data/analysis/ENm008/",CELL,".p1500_d3300.c840.a5.sis.node.dst",sep="")
	
}

ContactPair = read.table(ContactFile)[,1:2]

ListFile = paste(DataDir,"/List.txt",sep="")
FNWeight = read.table(ListFile)

# MaxNodes = 54
# NumRow = nrow(FNWeight)
# NumCol = MaxNodes*(MaxNodes+1)/2-MaxNodes - (MaxNodes-1)
# 
# DistIJMatrix = matrix(0,nrow=NumRow,ncol=NumCol)
# iCount = 1
# for (i in 1:(MaxNodes-2)){
# 	for (j in (i+2):MaxNodes){
# 		cat(i,j,"\n")
# 		flush.console()
# 		DistFile = paste(DataDir,"/",i,"_",j,".txt",sep="")
# 		DistIJMatrix[,iCount] = read.table(DistFile)[,1]
# 		iCount = iCount + 1
# 	}
# }


NumRow = nrow(FNWeight)
NumCol = nrow(ContactPair)

DistIJMatrix = matrix(0,nrow=NumRow,ncol=NumCol)
iCount = 1

for (irow in 1:nrow(ContactPair)){
	i = ContactPair[irow,1]
	j = ContactPair[irow,2]
	cat(i,j,"\n")
	flush.console()
	DistFile = paste(DataDir,"/",i,"_",j,".txt",sep="")
	DistIJMatrix[,iCount] = read.table(DistFile)[,1]
	iCount = iCount + 1
}

NumRecord = 300
NumRecord = NumRow
ptm <- proc.time()
RMSDMatrix = dist(DistIJMatrix[1:NumRecord,])/sqrt(NumCol)
proc.time() - ptm

GetClusterRatio <- function(DS, Weight){
	# DS = ds
	# Weight = FNWeight[1:NumRecord,2]
	
	Weight.Cluster = cbind(Weight,DS$cluster)
	TotalWeight = sum(Weight)

	ClusterIndexList = unique(sort(DS$cluster))
	TotalWeight.Cluster = matrix(,nrow=3,ncol=length(ClusterIndexList))
	iCount = 1
	for (i in ClusterIndexList){
		WeightSum = sum(Weight.Cluster[which(Weight.Cluster[,2]==i),1])
		WeightRatio = sprintf("%.2f%%",WeightSum*100/TotalWeight)
		TotalWeight.Cluster[,iCount] = c(i,WeightSum,WeightRatio)
		iCount = iCount + 1
	}
	Weight.DF = data.frame(cluster=ClusterIndexList,weight = as.numeric(TotalWeight.Cluster[2,]),ratio = as.character(TotalWeight.Cluster[3,]),stringsAsFactors = FALSE)
	return(Weight.DF)
}

PrintDBSCAN <- function(x, weight, ...){
	cat("dbscan Pts=", length(x$cluster), " MinPts=", x$MinPts, 
	    " eps=", x$eps, "\n", sep = "")
	if (is.null(x$isseed)) 
	    tab <- table(x$cluster)
	else {
	    tab <- table(c("seed", "border")[2 - x$isseed], cluster = x$cluster)
	    if (is.null(dim(tab))) {
	        tab <- cbind(tab)
	        colnames(tab) <- unique(x$cluster)
	    }
	    tab <- rbind(tab, total = colSums(tab))
	}
	
	print(tab, ...)
	weight.df = GetClusterRatio(x,weight)
	print(weight.df,...)
}

if (CELL == "GM12878"){
	EPSList = c(30:43)*10
}else{
	EPSList = c(20:29)*10
}

for (EPS in EPSList){
	ds = dbscan(RMSDMatrix,eps=EPS,MinPts=5,method="dist")
	OutFile = paste("/tmp/",CELL,".eps.",EPS,".cluster.txt",sep="")
	capture.output(file=OutFile,PrintDBSCAN(ds,weight=FNWeight[1:NumRecord,2]))
	write.table(
		file=OutFile,
		cbind(
			as.character(FNWeight[1:NumRecord,1]),
			FNWeight[1:NumRecord,2],
			ds$cluster,as.numeric(ds$isseed)),
		append=T,quote=F,row.names=F,col.name=F)	
}

