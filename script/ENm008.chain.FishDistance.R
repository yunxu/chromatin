#!/usr/bin/R
#------------------------------------------------------------------
# ENm008.chain.ClusterPvalue.R
#------------------------------------------------------------------
rm(list=ls())
# library(fpc)
#------------------------------------------------------------------

GM_cell = "GM12878"
K_cell = "K562"
Fish_Seg_1 = c(29779,86128)
Fish_Seg_2 = c(380160,433293)

#------------------------------------------------------------------
Mean_Pos_1 = mean(Fish_Seg_1)
Mean_Pos_2 = mean(Fish_Seg_2)
#------------------------------------------------------------------
All_Points_Pos = c(1,5693,15091,18344,29756,44231,50868,55911,64056,74448,88084,95256,100530,104687,109838,123474,131220,134334,147970,161606,167103,171769,185994,189074,203353,217802,225341,238977,247214,260850,274486,277942,289198,303867,310528,314538,327240,334889,348525,352385,360924,372670,380160,393796,407432,418222,421291,433293,445126,454365,468001,483412,496513,499411)
GetGeneNameAndRange <- function(){
  FN = "data/ENm008/circos/text.genes.txt"
  GeneRangeName = read.table(FN)[2:4]
  PosToName = NULL
  for (i in 1:length(All_Points_Pos)){
    pos = All_Points_Pos[i]
    str = ""
    for (j in 1:nrow(GeneRangeName)){
      if (  GeneRangeName[j,1] <= pos && GeneRangeName[j,2] >= pos ){
        if (str == ""){
          str = as.character(GeneRangeName[j,3])
        } else{
          str = paste(str, as.character(GeneRangeName[j,3]),sep=" / ")
        }
      }
    }
    # cat (i,"\t",str,"\n")
    PosToName = rbind(PosToName,c(i,str))
  }
  return(PosToName)
}

PosToName = GetGeneNameAndRange()
xlabel = apply(PosToName,1,function(x){paste(x[2],x[1])})

#----------------------------------------------------------------------------
GetClosestEquPointIndFromExpPointPos <- function(ExpPointPos){
  ExpPointInd = NULL
  for (i in 1:length(ExpPointPos)){
    Points_Diff = abs(All_Points_Pos - ExpPointPos[i])
    ExpPointInd = c(ExpPointInd,which(Points_Diff == min(Points_Diff)))
  }
  return (ExpPointInd)
}

#----------------------------------------------------------------------------
Node_1 = GetClosestEquPointIndFromExpPointPos(Mean_Pos_1)
Node_2 = GetClosestEquPointIndFromExpPointPos(Mean_Pos_2)

DataPath = "/Users/yunxu/workspace/projects/chromatin/result/analysis/ENm008/chain/weight.dist/par.1_1_1_1/"
FN_GM = paste(DataPath,GM_cell,".p1500_d3300/",Node_1,"_",Node_2,".txt",sep="")
Dist_Weight_GM = read.table(FN_GM)
Dist_GM = rep(Dist_Weight_GM[,1],Dist_Weight_GM[,2])
Dist_GM = Dist_GM / 10 # A to nm
FN_K = paste(DataPath,K_cell,".p1500_d3300/",Node_1,"_",Node_2,".txt",sep="")
Dist_Weight_K = read.table(FN_K)
Dist_K = rep(Dist_Weight_K[,1],Dist_Weight_K[,2])
Dist_K = Dist_K / 10 # A to nm

#----------------------------------------------------------------------------
library(reshape2)
library(ggplot2)
library(RColorBrewer)

#GM.Fish.Exp = data.frame("cell"= "GM12878", "mean" = 320, "sem"= 20)
#K.Fish.Exp = data.frame("cell" = "K562", "mean" = 390, "sem" = 20)
GM.Fish.Exp = data.frame("cell"= "GM12878", "mean" = 320, "sem"= NA)
K.Fish.Exp = data.frame("cell" = "K562", "mean" = 390, "sem" = NA)
Fish.Exp = rbind(GM.Fish.Exp, K.Fish.Exp)
Fish.Exp$type = "Exp"

#GM.Fish.Prd = data.frame("cell" = "GM12878", "mean"= mean(Dist_GM), "sem" = sd(Dist_GM/length(Dist_GM))) 
#K.Fish.Prd = data.frame("cell" = "K562", "mean"= mean(Dist_K), "sem" = sd(Dist_K)/length(Dist_K))

GM.Fish.Prd = data.frame("cell" = "GM12878", "mean"= mean(Dist_GM), "sem" = sd(Dist_GM)) 
K.Fish.Prd = data.frame("cell" = "K562", "mean"= mean(Dist_K), "sem" = sd(Dist_K))
Fish.Prd = rbind(GM.Fish.Prd, K.Fish.Prd)
Fish.Prd$type = "Prd"

Fish.Data = rbind(Fish.Exp, Fish.Prd)
MyPal = brewer.pal(9,"Set1")
dodge <- position_dodge(width=.8)
limits <- aes(ymin = mean - sem, ymax = mean + sem)

gp <- ggplot(Fish.Data, aes(x=type,y=mean,fill=cell, width=.5 ))
gp <- gp +
   theme_bw() +
  theme(text = element_text(size=20),
    legend.justification=c(0,1),
    legend.direction = "horizontal",
    legend.position=c(0,1) ) +
  labs(x= paste(""),y = "Distance (nm)") +
  scale_x_discrete(label=c("FISH","Prediction")) +
  coord_cartesian(ylim=c(0,500)) +
  scale_fill_manual("",
    values = c( "GM12878" = MyPal[1], "K562" = MyPal[2]))
 
gp = gp + geom_bar(width=0.7, stat="identity",position=dodge)  + geom_errorbar(limits, position=dodge, width=0.25)

postscript(file="/tmp/fish.eps", width=6, height=4)
  gp
dev.off()

q()

#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------

GM.df = melt(GM.DistMat)
GM.df$cell = "GM"
K.df = melt(K.DistMat)
K.df$cell = "K"
df = rbind(GM.df[,2:ncol(GM.df)],K.df[,2:ncol(K.df)])
colnames(df) = c("NodeNode","distance","cell")
# colnames(K.df) = c("NodeNode","distance","cell")
ggp <- ggplot(df,aes(NodeNode,distance)) + 
	theme_bw() +
	theme(axis.text.x=element_text(angle=90, hjust=0)) +
	geom_boxplot(aes(fill=cell),outlier.shape = NA) + 
	scale_fill_manual(name="cell",values=c(brewer.pal(9,"Set1")	)) +	
	# scale_colour_brewer(palette="Set1") +
	labs(x= paste("Node - Node"),y = "distance (nm)") +
  scale_y_continuous(breaks=seq(0,900,by=100)) 
ggp

#----------------------------------------------------------------------------



#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
Node_1_vec = GetClosestEquPointIndFromExpPointPos(Fish_Seg_1)
Node_2_vec = GetClosestEquPointIndFromExpPointPos(Fish_Seg_2)

Grp_1_List = seq(Node_1_vec[1],Node_1_vec[2])
Grp_2_List = seq(Node_2_vec[1],Node_2_vec[2])

Pair_List = expand.grid(Grp_1_List,Grp_2_List)

#----------------------------------------------------------------------------
SampleFile_GM = paste(DataPath,GM_cell,".p1500_d3300/",Pair_List[1,1],"_",Pair_List[1,2],".txt",sep="") 
SampleData_GM = read.table(SampleFile_GM)
NRow_GM = nrow(SampleData_GM)
Weight_GM = SampleData_GM[,2]
NumPair = nrow(Pair_List)

DistMat_GM = matrix(0,nrow=NRow_GM,ncol=NumPair)
for (i in 1:nrow(Pair_List)){
  Node_1 = Pair_List[i,1]
  Node_2 = Pair_List[i,2]
  FN_GM = paste(DataPath,GM_cell,".p1500_d3300/",Node_1,"_",Node_2,".txt",sep="") 
  Dist_Weight_GM = read.table(FN_GM)
  DistMat_GM[,i] = Dist_Weight_GM[,1]/10 
}
summary(apply(DistMat_GM,1,mean))



SampleFile_K = paste(DataPath,K_cell,".p1500_d3300/",Pair_List[1,1],"_",Pair_List[1,2],".txt",sep="") 
SampleData_K = read.table(SampleFile_K)
NRow_K = nrow(SampleData_K)
Weight_K = SampleData_K[,2]
NumPair = nrow(Pair_List)

DistMat_K = matrix(0,nrow=NRow_K,ncol=NumPair)
for (i in 1:nrow(Pair_List)){
  Node_1 = Pair_List[i,1]
  Node_2 = Pair_List[i,2]
  FN_K = paste(DataPath,K_cell,".p1500_d3300/",Node_1,"_",Node_2,".txt",sep="") 
  Dist_Weight_K = read.table(FN_K)
  DistMat_K[,i] = Dist_Weight_K[,1]/10 
}
summary(apply(DistMat_K,1,mean))
# > summary(apply(DistMat_GM,1,mean))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   99.77  125.10  141.30  145.10  160.50  254.70 
# > summary(apply(DistMat_K,1,mean))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   133.3   248.3   308.9   313.1   372.6   645.0 
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
GetWeight <- function(ListFN){
	NumFiles = length(ListFN)
	LogWeight = rep(0, NumFiles)
	for (i in 1:NumFiles){
		if (i %% 100 == 0){cat (i,"/",NumFiles,"\n"); flush.console()}
		FN = ListFN[i]
		LogWeight[i] = read.table(FN,comment.char="",nrow=1)[,3]
	}
	return(LogWeight)
}
GetTwoGrpDist <- function(ListFN,Grp1,Grp2){
	NumFiles = length(ListFN)
	Dist = rep(0,NumFiles)
	for (i in 1:NumFiles){
		if (i %% 100 == 0){cat (i,"/",NumFiles,"\n"); flush.console()}
		FN = ListFN[i]
		Pts = read.table(FN)
		Pts = Pts/10
		AvgPts1 = Pts[Grp1,1:3]
		AvgPts2 = Pts[Grp2,1:3]
		Dist[i] = sqrt(sum((colMeans(AvgPts1)-colMeans(AvgPts2))^2))
	}
	return(Dist)
}
GetListGrpDist <- function(ListFN,CountVec,GrpList){
	# ListFN = GM.ListFN
	# CountVec = GM.Count
	NumFiles = length(ListFN)
	# NumFiles = 100
	# Dist = rep(0,NumFiles)
	
	DistMat = matrix(0,nrow=NumFiles,ncol=nrow(GrpList))
	for (i in 1:NumFiles){
		if (i %% 100 == 0){cat (i,"/",NumFiles,"\n"); flush.console()}
		FN = ListFN[i]
		Pts = read.table(FN)[,1:3]
		Pts = Pts/10
		PtsDiff = Pts[GrpList[,1],] - Pts[GrpList[,2],]
		DistMat[i,] = sqrt(rowSums(PtsDiff * PtsDiff))
	}
	RepDistMat = apply(DistMat,2,function(x){rep(x,CountVec[1:NumFiles])})
	colnames(RepDistMat) = apply(GrpList,1,function(x){paste(x,collapse="-")})
	return(RepDistMat)
}

#----------------------------------------------------------------------------
PtsPath_GM = paste("/Users/yunxu/workspace/projects/chromatin/result/ENm008/chain/par.1_1_1_1/",GM_cell,".p1500_d3300",sep="")
FileList_GM = dir(PtsPath_GM,recursive=T,pattern=".pts",full.names=T)

PtsPath_K = paste("/Users/yunxu/workspace/projects/chromatin/result/ENm008/chain/par.1_1_1_1/",K_cell,".p1500_d3300",sep="")
FileList_K = dir(PtsPath_K,recursive=T,pattern=".pts",full.names=T)

GM.LogWeight = GetWeight(FileList_GM)
K.LogWeight = GetWeight(FileList_K)

GM.Weight = exp(GM.LogWeight-max(GM.LogWeight))
K.Weight = exp(K.LogWeight - max(K.LogWeight))

NumSamples = 1000
GM.Count = round(NumSamples*GM.Weight)
K.Count = round(NumSamples*K.Weight)

GM.Ind = which(GM.Count!=0)
K.Ind = which(K.Count != 0)


GM.Count = GM.Count[GM.Ind]
K.Count = K.Count[K.Ind]

GM.ListFN = FileList_GM[GM.Ind]
K.ListFN = FileList_K[K.Ind]

GM.DistMat = GetListGrpDist(GM.ListFN,GM.Count,Pair_List)
K.DistMat = GetListGrpDist(K.ListFN,K.Count,Pair_List)
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
#----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

GM.ChainListFN.File = "/tmp/GM12878.chain.ls.txt"
K.ChainListFN.File = "/tmp/K562.chain.ls.txt"

GM.ListFN = as.character(read.table(GM.ChainListFN.File)[,1])
K.ListFN = as.character(read.table(K.ChainListFN.File)[,1])

GetWeight <- function(ListFN){
	NumFiles = length(ListFN)
	LogWeight = rep(0, NumFiles)
	for (i in 1:NumFiles){
		if (i %% 100 == 0){cat (i,"/",NumFiles,"\n"); flush.console()}
		FN = ListFN[i]
		LogWeight[i] = read.table(FN,comment.char="",nrow=1)[,3]
	}
	return(LogWeight)
}

GetTwoGrpDist <- function(ListFN,Grp1,Grp2){
	NumFiles = length(ListFN)
	Dist = rep(0,NumFiles)
	for (i in 1:NumFiles){
		if (i %% 100 == 0){cat (i,"/",NumFiles,"\n"); flush.console()}
		FN = ListFN[i]
		Pts = read.table(FN)
		Pts = Pts/10
		AvgPts1 = Pts[Grp1,1:3]
		AvgPts2 = Pts[Grp2,1:3]
		Dist[i] = sqrt(sum((colMeans(AvgPts1)-colMeans(AvgPts2))^2))
	}
	return(Dist)
}

GM.LogWeight = GetWeight(GM.ListFN)
K.LogWeight = GetWeight(K.ListFN)

GM.Weight = exp(GM.LogWeight-max(GM.LogWeight))
K.Weight = exp(K.LogWeight - max(K.LogWeight))

NumSamples = 1000
GM.Count = round(NumSamples*GM.Weight)
K.Count = round(NumSamples*K.Weight)

GM.Ind = which(GM.Count!=0)
K.Ind = which(K.Count != 0)

GM.Count = GM.Count[GM.Ind]
K.Count = K.Count[K.Ind]
GM.ListFN = GM.ListFN[GM.Ind]
K.ListFN = K.ListFN[K.Ind]

Grp1 = c(4:12)
# Grp1 = c(4:5)
Grp1 = c(5)
Grp2 = c(42:48)
# Grp2 = c(42:43)
Grp2 = c(47)


GrpList = expand.grid(Grp1,Grp2)

GetListGrpDist <- function(ListFN,CountVec,GrpList){
	# ListFN = GM.ListFN
	# CountVec = GM.Count
	NumFiles = length(ListFN)
	# NumFiles = 100
	# Dist = rep(0,NumFiles)
	
	DistMat = matrix(0,nrow=NumFiles,ncol=nrow(GrpList))
	for (i in 1:NumFiles){
		if (i %% 100 == 0){cat (i,"/",NumFiles,"\n"); flush.console()}
		FN = ListFN[i]
		Pts = read.table(FN)[,1:3]
		Pts = Pts/10
		PtsDiff = Pts[GrpList[,1],] - Pts[GrpList[,2],]
		DistMat[i,] = sqrt(rowSums(PtsDiff * PtsDiff))
	}
	RepDistMat = apply(DistMat,2,function(x){rep(x,CountVec[1:NumFiles])})
	colnames(RepDistMat) = apply(GrpList,1,function(x){paste(x,collapse="-")})
	return(RepDistMat)
}

GM.DistMat = GetListGrpDist(GM.ListFN,GM.Count,GrpList)
K.DistMat = GetListGrpDist(K.ListFN,K.Count,GrpList)

library(reshape2)
library(ggplot2)
library(RColorBrewer)
GM.df = melt(GM.DistMat)
GM.df$cell = "GM"
K.df = melt(K.DistMat)
K.df$cell = "K"
df = rbind(GM.df[,2:ncol(GM.df)],K.df[,2:ncol(K.df)])
colnames(df) = c("NodeNode","distance","cell")
# colnames(K.df) = c("NodeNode","distance","cell")
ggp <- ggplot(df,aes(NodeNode,distance)) + 
	theme_bw() +
	theme(axis.text.x=element_text(angle=90, hjust=0)) +
	geom_boxplot(aes(fill=cell),outlier.shape = NA) + 
	scale_fill_manual(name="cell",values=c(brewer.pal(9,"Set1")	)) +	
	# scale_colour_brewer(palette="Set1") +
	labs(x= paste("Node - Node"),y = "distance (nm)") +
  scale_y_continuous(breaks=seq(0,900,by=100)) 
ggp



ggp1 <- ggplot(df,aes(x=NodeNode,y=distance)) +
	theme_bw() +
	geom_boxplot(aes(fill=cell)) + 
	scale_fill_manual(name="cell",values=c(brewer.pal(9,"Set1")	)) +	
	labs(x= paste("Node - Node"),y = "distance (nm)") +
	labs(title=paste("Node",paste(GrpList,collapse="-"))) +
  scale_y_continuous(breaks=seq(0,900,by=100)) 
ggp2 <- ggplot(df,aes(x=distance)) +
	theme_bw() +
	geom_histogram(aes(fill=cell),binwidth=10) +
	scale_fill_manual(name="cell",values=c(brewer.pal(9,"Set1")	)) +	
	labs(title=paste("Node",paste(GrpList,collapse="-"))) +
	scale_x_continuous(breaks=seq(0,900,by=100)) +
	facet_grid(cell~.)
library(gridExtra)
grid.arrange(ggp1,ggp2,nrow=1)




# GM.Mat = matrix(0,nrow=sum(GM.Count),ncol=nrow(Grp),
# 	dimnames=list(NULL,apply(Grp,1,function(x){paste(x,collapse="-")})))
# K.Mat = matrix(0,nrow=sum(K.Count),ncol=nrow(Grp),
# 	dimnames=list(NULL,apply(Grp,1,function(x){paste(x,collapse="-")})))
# for ( i in 1:nrow(Grp)){
# 	Grp1 = Grp[i,1]
# 	Grp2 = Grp[i,2]
# 	GM.Grp.Dist = GetTwoGrpDist(GM.ListFN ,Grp1,Grp2)
# 	K.Grp.Dist = GetTwoGrpDist(K.ListFN ,Grp1,Grp2)
# 	GM.Dist.Samples = rep(GM.Grp.Dist,GM.Count)
# 	K.Dist.Samples = rep(K.Grp.Dist,K.Count)
# 	NodeNodeStr = paste(Grp[i,],collapse="-")
# 	GM.Mat[,i] = GM.Dist.Samples
# 	K.Mat[,i] = K.Dist.Samples
# }
# 
# GM.df <- data.frame("dist" = GM.Dist.Samples,"cell" = as.factor("GM"),"node-node"=NodeNodeStr)
# K.df <- data.frame("dist" = K.Dist.Samples,"cell"= as.factor("K"),"node-node"=NodeNodeStr)
# 
#------------------------------------------------------------------

Grp1 = c(4:12)
Grp2 = c(42:48)
Grp1 = Node_1_vec
Grp2 = Node_2_vec
GM.Grp.Dist = GetTwoGrpDist(GM.ListFN,Grp1,Grp2)
K.Grp.Dist = GetTwoGrpDist(K.ListFN,Grp1,Grp2)
# 
GM.Dist.Samples = rep(GM.Grp.Dist,GM.Count)
K.Dist.Samples = rep(K.Grp.Dist,K.Count)
GM.df <- data.frame("distance" = GM.Dist.Samples,"cell" = as.factor("GM"))
K.df <- data.frame("distance" = K.Dist.Samples,"cell"= as.factor("K"))
# 
df <- rbind(GM.df,K.df)

ggp1 <- ggplot(df,aes(x=cell,y=distance)) +
	theme_bw() +
	geom_boxplot(aes(fill=cell)) + 
	scale_fill_manual(name="cell",values=c(brewer.pal(9,"Set1")	)) +	
	labs(x= paste("Nodes - Nodes"),y = "distance (nm)") +
	labs(title= paste(paste(c(range(Grp1)),collapse=":"),paste(c(range(Grp2)),collapse=":"),sep="-") ) +
  scale_y_continuous(breaks=seq(0,900,by=100)) 
ggp2 <- ggplot(df,aes(x=distance)) +
	theme_bw() +
	geom_histogram(aes(fill=cell),binwidth=10) +
	scale_fill_manual(name="cell",values=c(brewer.pal(9,"Set1")	)) +	
	labs(x="distance (nm)") +
	scale_x_continuous(breaks=seq(0,900,by=100)) +
	facet_grid(cell~.)
library(gridExtra)
grid.arrange(ggp1,ggp2,nrow=1)

# ggp <- ggplot(df,aes(x=cell,y=dist)) + geom_boxplot() + theme_bw()
# 
# boxplot(list(GM.Dist.Samples,K.Dist.Samples))
