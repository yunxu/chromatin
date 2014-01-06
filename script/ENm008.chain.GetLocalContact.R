#!/usr/bin/R

rm(list=ls())
#------------------------------------------------------------------
library(ggplot2)
library(scales)
#------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
Type="r270.ovl"
method="alpha"
if (Sys.info()["sysname"] == "Darwin"){
	CELL = "K562"
	CELL = "GM12878"
	ConDir=paste("/Users/yunxu/workspace/projects/chromatin/result/analysis/ENm008/chain/",method,"/par.1_1_1_1/",CELL,".p1500_d3300",sep="")
}else{
	CELL = args[1]
	ConDir=paste("/dump/yunxu/working/projects.new/chromatin/result/analysis/ENm008/chain/",method,"/par.1_1_1_1/",CELL,".p1500_d3300",sep="")
}
print(ConDir)
#------------------------------------------------------------------
# NumMax = 1000
OutDir = "tmp/"
OutDir = "local/"
#------------------------------------------------------------------
#------------------------------------------------------------------
SegLenFN="data/analysis/ENm008/ENm008.p1500.equ.seg.len.txt"
NumNodes=nrow(read.table(SegLenFN))+1
#------------------------------------------------------------------

GetWeightFile <- function(ConDir){
	
	ConFiles   = list.files(path=ConDir,pattern=Type,full.names=T)
	NumFiles = length(ConFiles)

	BaseName = sub("[.][^.]*$", "", ConFiles)
	OverlapFiles = paste(BaseName,".05.a.overlap",sep="")
# NumFiles = 1000
	Weight = rep(0,NumFiles)
	ConList = list()
	for (i in 1:NumFiles){
		Weight[i] = read.table(ConFiles[i],comment.char="",nrows=1)[,3]
		ConList[[i]] = read.table(OverlapFiles[i])[,1:2]
	}
	return(list("weight"=Weight,"conlist"=ConList))
}
#------------------------------------------------------------------
GetSameConnection <- function(Con,ConWeight,Node1,Node2){
	# Node1 <-> Node2
	#     \    /
	#       k
	# Con = ConList[[1]]
	# Node1 = 14
	# Node2 = 27
	HaveConnection = which(Con[,1] == Node1 & Con[,2] == Node2)
	LinkWeight = 0
	SameConnections = integer(0)
	if (length(HaveConnection)){
		LinkWeight = ConWeight
		ConnectWithNode1 = Con[which(Con[,1]==Node1),2]
		ConnectWithNode2 = Con[which(Con[,1]==Node2),2]
		ConnectWithNode1and2 = intersect(ConnectWithNode1, ConnectWithNode2)
		# remove step=1 neighbour
		NeighbourStep1 = which(abs(ConnectWithNode1and2-Node1)==1 | abs(ConnectWithNode1and2-Node2)==1)
		SameConnections = setdiff(ConnectWithNode1and2,ConnectWithNode1and2[NeighbourStep1])
	}
	return(list("linkweight"=LinkWeight, "con" = SameConnections))
}
#------------------------------------------------------------------

GetNodeInd_Node1Node2Connect <- function(Con,ConWeight,Node1,Node2,NodeIndSel){
	# K <-> Node1 <-> Node2 or 
	# Node1 <-> Node2 <-> K
	
# Under condition Node1-Node2	
# Con = ConList[[1]]
# Node1 = 14
# Node2 = 27
# NodeIndSel = 14
	HaveConnection = which(Con[,1] == Node1 & Con[,2] == Node2)
	LinkWeight = 0
	NodeInd = integer(0)
	if (length(HaveConnection)){
		ConnectWithNode1 = Con[which(Con[,1]==Node1),2]
		ConnectWithNode2 = Con[which(Con[,1]==Node2),2]
		
		if (NodeIndSel == Node1){
			NodeInd = setdiff(ConnectWithNode1, ConnectWithNode2)
			NodeInd = NodeInd[which(NodeInd != Node2)]
			NeighbourStep1 = which(abs(NodeInd-Node1)==1)
		}else if(NodeIndSel == Node2){
			NodeInd = setdiff(ConnectWithNode2, ConnectWithNode1)
			NodeInd = NodeInd[which(NodeInd != Node1)]
			NeighbourStep1 = which(abs(NodeInd-Node2)==1)
		}
		# remove step=1 neighbour
		# NeighbourStep1 = which(abs(NodeInd-Node1)==1 | abs(NodeInd-Node2)==1)
		NodeInd = setdiff(NodeInd,NodeInd[NeighbourStep1])
		if (length(NodeInd)){
			LinkWeight = ConWeight
		}
	}
	return(list("linkweight"=LinkWeight,"con"=NodeInd))
}
#------------------------------------------------------------------

GetNodeInd_Node1Node2NotConnect <- function(Con,ConWeight,Node1,Node2){ 
	# Node1<->K<->Node2
	# Con = ConList[[1]]
	# Node1 = 14
	# Node2 = 40
	HaveConnection = which(Con[,1] == Node1 & Con[,2] == Node2)
	LinkWeight = 0
	SameConnections = integer(0)
	if (length(HaveConnection) == 0){ # no connection
		ConnectWithNode1 = Con[which(Con[,1]==Node1),2]
		ConnectWithNode2 = Con[which(Con[,1]==Node2),2]
		ConnectWithNode1and2 = intersect(ConnectWithNode1, ConnectWithNode2)
		# remove step=1 neighbour
		NeighbourStep1 = which(abs(ConnectWithNode1and2-Node1)==1 | abs(ConnectWithNode1and2-Node2)==1)
		SameConnections = setdiff(ConnectWithNode1and2,ConnectWithNode1and2[NeighbourStep1])
		if(length(SameConnections)){
			LinkWeight = ConWeight
		}
	}
	return (list("linkweight"=LinkWeight,"con"=SameConnections))
}
#------------------------------------------------------------------

GetWeightCon <- function(ConList,LinkWeight.Vec){
	# ConList = SameConnectionList
	ConInd = unique(sort(unlist(ConList)))
	NumCon = length(ConInd)
	WeightMat = matrix(0,nrow=NumCon,ncol=2)
	if(length(ConInd)){
		for (i in 1:NumCon){
			Con = unlist(lapply(ConList,function(x){length(which(x==ConInd[i]))}))
			NormalizedWeight = sum(LinkWeight.Vec * Con) / sum(LinkWeight.Vec)
			WeightMat[i,] = c(ConInd[i],NormalizedWeight)
		}
	}
	return (WeightMat)
}
#------------------------------------------------------------------


GetResult <- function(Node1,Node2){
	# Node1 = 3
	# Node2 = 5
	SameConnectionList = list()
	SameConnectionWeight = rep(0,NumFiles)
	Node1ConnectionList = list()
	Node1ConnectionWeight = rep(0,NumFiles)
	Node2ConnectionList = list()
	Node2ConnectionWeight = rep(0,NumFiles)
	IndirectConnectionList = list()
	IndirectConnectionWeight = rep(0,NumFiles)

	for (i in 1:NumFiles){
		# Node1 = 14
		# Node2 = 20
		# # Node2 = 27
		listtmp = GetSameConnection(ConList[[i]],Weight[i],Node1,Node2)
		SameConnectionList[[i]] = listtmp$con
		SameConnectionWeight[i] = listtmp$linkweight
		listtmp = GetNodeInd_Node1Node2Connect(ConList[[i]],Weight[i],Node1,Node2,Node1)
		Node1ConnectionList[[i]] = listtmp$con
		Node1ConnectionWeight[i] = listtmp$linkweight
		listtmp = GetNodeInd_Node1Node2Connect(ConList[[i]],Weight[i],Node1,Node2,Node2)
		Node2ConnectionList[[i]] = listtmp$con
		Node2ConnectionWeight[i] = listtmp$linkweight
		listtmp = GetNodeInd_Node1Node2NotConnect(ConList[[i]],Weight[i],Node1,Node2)
		IndirectConnectionList[[i]] = listtmp$con
		IndirectConnectionWeight[i] = listtmp$linkweight
	}

	SameConnectionInd = unique(sort(unlist(SameConnectionList)))
	Node1ConnectionInd = unique(sort(unlist(Node1ConnectionList)))
	Node2ConnectionInd = unique(sort(unlist(Node2ConnectionList)))
	IndirectConnectionInd = unique(sort(unlist(IndirectConnectionList)))

	

	SameConnectionWeightMat = GetWeightCon(SameConnectionList,SameConnectionWeight)
	Node1ConnectionWeightMat = GetWeightCon(Node1ConnectionList,Node1ConnectionWeight)
	Node2ConnectionWeightMat = GetWeightCon(Node2ConnectionList,Node2ConnectionWeight)
	IndirectConnectionWeightMat = GetWeightCon(IndirectConnectionList,IndirectConnectionWeight)

	if (nrow(SameConnectionWeightMat)==0){
    stop(paste(Node1,"-",Node2,"has no connection\n"))
	}
	SameDF = data.frame(NodeIndex=SameConnectionWeightMat[,1],value=SameConnectionWeightMat[,2],variable="Triplet")
	Node1DF = data.frame(NodeIndex=Node1ConnectionWeightMat[,1],value=Node1ConnectionWeightMat[,2],variable=paste("Connect-",Node1,sep=""))
	Node2DF = data.frame(NodeIndex=Node2ConnectionWeightMat[,1],value=Node2ConnectionWeightMat[,2],variable=paste("Connect-",Node2,sep=""))
	IndirectDF =data.frame(NodeIndex=IndirectConnectionWeightMat[,1],value=IndirectConnectionWeightMat[,2],variable="Intercept")
	WeightDF= rbind(SameDF,Node1DF,Node2DF,IndirectDF)

#------------------------------------------------------------------
	p <- ggplot(data=WeightDF,aes(x=NodeIndex,y=value,color=variable)) +
		theme_bw() +
		geom_line() + 
		geom_point() +
		theme(legend.position = "none") +
		labs(title=CELL, x= paste("Node Index"),y = paste("Node (",Node1 , "-", Node2,") Contact Propensity")) + 
		scale_colour_brewer(palette="Set1") +
		scale_x_continuous(breaks=seq(5,NumNodes,by=5)) +
		scale_y_continuous(breaks=seq(0,0.6,by=0.2))
	vline.data1 <- data.frame(variable = levels(WeightDF$variable), vl=c(Node1,Node1,Node2,Node1)) 
	vline.data2 <- data.frame(variable = levels(WeightDF$variable), vl=c(Node2,Node1,Node2,Node2)) 
	# dev.new(width=3,height=6)
	p = p + facet_grid(variable ~.) + 
		geom_vline(aes(xintercept = vl), vline.data1,linetype="dashed",size=.5,colour="grey40") +
		geom_vline(aes(xintercept = vl), vline.data2,linetype="dashed",size=.5,colour="grey40") +
		coord_cartesian(ylim=c(0,0.8))

	(FN <- paste(OutDir,"/",CELL,".p1500_d3300.c840.",Type,".local.",Node1,"-",Node2,".pdf",sep=""))
	pdf(file=FN,width=3,height=5)
	print(p)
	dev.off()
#------------------------------------------------------------------
	WM = matrix(0,nrow=NumNodes,ncol=5)
	WM[,1] = 1:NumNodes
	WM[SameConnectionWeightMat[,1],2] = SameConnectionWeightMat[,2]
	WM[Node1ConnectionWeightMat[,1],3] = Node1ConnectionWeightMat[,2]
	WM[Node2ConnectionWeightMat[,1],4] = Node2ConnectionWeightMat[,2]
	WM[IndirectConnectionWeightMat[,1],5] = IndirectConnectionWeightMat[,2]
	colnames(WM) = c("NodeIndex", "Triplet", "K-Node1", "Node2-K", "Intercept")
	
	(FN <- paste(OutDir,"/",CELL,".p1500_d3300.c840.",Type,".local.",Node1,"-",Node2,".txt",sep=""))
	write.table( file=FN,
		paste("#",paste(colnames(WM),collapse="\t")), 
		quote=F,row.names=F,col.names=F)
	write.table( file =FN,
		apply(WM,1,function(x){sprintf("%g\t%.3f\t%.3f\t%.3f\t%.3f",x[1],x[2],x[3],x[4],x[5])}),
		quote=F,row.names=F,col.names=F,append=T)
	
}

WeightData <- GetWeightFile(ConDir)
Weight <- WeightData$weight
ConList <- WeightData$conlist
NumFiles <- length(Weight)

InterestingPair = NULL
for(i in 3:(NumNodes-3)){
	for (j in (i+3):NumNodes){
		InterestingPair = rbind(InterestingPair, c(i,j))
	}
}

InterestingPair = NULL
InterestingNodes = c(3,7,9,10,13,17,18,21,34,35,36,43,52)
#InterestingNodes = c(9,10,13,17,18,21,34,35,36,43,52)
InterestingPair = t(combn(InterestingNodes,2))
InterestingPair = InterestingPair[which(InterestingPair[,2]-InterestingPair[,1] !=1),]
#print(InterestingPair)
#InterestingPair = NULL
#InterestingPair = rbind(InterestingPair,c(3,13))
#InterestingPair = rbind(InterestingPair,c(7,13))
#InterestingPair = rbind(InterestingPair,c(9,13))
#InterestingPair = rbind(InterestingPair,c(10,13))
#InterestingPair = rbind(InterestingPair,c(13,17))
#InterestingPair = rbind(InterestingPair,c(13,18))
#InterestingPair = rbind(InterestingPair,c(13,21))
#InterestingPair = rbind(InterestingPair,c(13,34))
#InterestingPair = rbind(InterestingPair,c(13,35))
#InterestingPair = rbind(InterestingPair,c(13,36))
#InterestingPair = NULL
#InterestingPair = rbind(InterestingPair,c(3,35))

for (i in 1:nrow(InterestingPair)){
	cat("(", paste(InterestingPair[i,],collapse=",") , ")\t", i,"/", nrow(InterestingPair),"\n")
	flush.console()
	GetResult(InterestingPair[i,1],InterestingPair[i,2])
}


