#!/usr/bin/R
###################################################################
# Get contact index according to segment (starting and end)
###################################################################

rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)

if (Sys.info()["sysname"] == "Darwin"){
	SampleSize    = "p1500_d3300"
}else{
	SampleSize    = args[1]
}
SampleSize    = "p1500_d3300"
#------------------------------------------------------------------
Type=".con"
RndConDir="../result/analysis/ENm008/ensemble/cut.all"
NewConDir_1 = "../result/analysis/ENm008/chain/cut.all/par.1_1_1_1/GM12878.p1500_d3300"
NewConDir_2 = "../result/analysis/ENm008/chain/cut.all/par.1_1_1_1/K562.p1500_d3300"
#------------------------------------------------------------------
#Type=".r150.ovl"
#Type=".r270.ovl"
#RndConDir="../result/analysis/ENm008/ensemble/alpha"
#NewConDir_1 = "../result/analysis/ENm008/chain/alpha/par.1_1_1_1/GM12878.p1500_d3300"
#NewConDir_2 = "../result/analysis/ENm008/chain/alpha/par.1_1_1_1/K562.p1500_d3300"
#------------------------------------------------------------------
# Parameter
#------------------------------------------------------------------
ReferenceCountFile = paste("data/analysis/ENm008/",SampleSize,".comb.c.txt",sep="")
cat(ReferenceCountFile, "\n") 
DestDir = "tmp/"
#------------------------------------------------------------------
# Read reference data, and merge to experimental segment
#------------------------------------------------------------------
RF.Vec  = read.table(ReferenceCountFile)
Con.Vec = RF.Vec[,1:2]

# GM12878.p1500_d3300.c840.a5.sis.node.dst
DstFile_1 = "data/analysis/ENm008/GM12878.p1500_d3300.c840.a5.sis.node.dst"
Con.Vec_1 = read.table(DstFile_1)[,1:2]
DstFile_2 = "data/analysis/ENm008/K562.p1500_d3300.c840.a5.sis.node.dst"
Con.Vec_2 = read.table(DstFile_2)[,1:2]

# # replace zero with minimal
# Contact.RF.Count.Vec = RF.Vec[,3:ncol(RF.Vec)]
# Contact.RF.Count.Vec =	apply(Contact.RF.Count.Vec,2,function(x)
# 			{ matchind = which(x==0)
# 				if(length(matchind)){x[matchind]= min(x[which(x!=0)])}
# 				return(x) })
# RF.Vec = cbind(Con.Vec, Contact.RF.Count.Vec)
# RF.Vec.Scale = cbind(Con.Vec, apply(Contact.RF.Count.Vec,2,function(x){x/sum(x)}))

# #------------------------------------------------------------------
GetWeight <- function(FNList){
	NumFiles = length(FNList)
	Val = rep (0,NumFiles)
	for (i in 1:NumFiles){
		FN = FNList[i]
		Val[i] = read.table(FN,comment.char="",nrows=1)[,3]
	}
	return(Val)
}

GetConFromFiles <- function(FNList){
	# FNList = NewConFiles[1,drop=F]
	NumFiles = length(FNList)
	ValMat = matrix(0, nrow=NumPairs,ncol=NumFiles)
	for (i in 1:NumFiles){
		# i = 1
		FN = FNList[i]
		if (i %% 100 == 0){
			cat(FN,"/",i,"/",NumFiles,"\n");flush.console();
		}
		Data = read.table(FN)
		ValMat[,i] = apply(AllCon.Vec,1,function(x)
			{
				length(which(x[1] == Data[,1] & x[2] == Data[,2]))	
			}	
		)
	}
	return(ValMat)
}

GenerateAllConVec <- function(NumNodes){
	AllCon.Vec = matrix(0,nrow=NumPairs,ncol=2)
	k = 1
	for (i in 3:NumNodes){
		for (j in 1:(i-2)){
			AllCon.Vec[k,] = c(j,i)
			k = k+1
		}
	}
	return(AllCon.Vec)
}

GetPIJFromMatWeight <- function(ConMat, ConWeight){
	# ConMat = NewConMat
	# ConWeight = NewConWeight
	SumWeightConMat = t(t(ConMat)*ConWeight)
	Sum_K_ConIJ = apply(SumWeightConMat,1,sum)
	Sum_KIJ_ConIJ = sum(SumWeightConMat)
	P_IJ = Sum_K_ConIJ/ Sum_KIJ_ConIJ
	return(P_IJ)
}

GetPIJFromDir <- function(ConDir){
	#ConFiles   = list.files(path=ConDir,pattern=".con",full.names=T)
	ConFiles   = list.files(path=ConDir,pattern=Type,full.names=T)
#	ConFiles   = ConFiles[1:10]
	ConWeight  = GetWeight(ConFiles)
	ConMat     = GetConFromFiles(ConFiles)
	ConPIJ     = GetPIJFromMatWeight(ConMat, ConWeight)
	return (list("weight" = ConWeight,"pij" = ConPIJ,"mat"=ConMat))
}

#------------------------------------------------------------------
# Get weight pij and matrix
#------------------------------------------------------------------
SegLenFN="data/analysis/ENm008/ENm008.p1500.equ.seg.len.txt"
NumNodes=nrow(read.table(SegLenFN))+1

#RndConFiles=list.files(path=RndConDir,pattern=".con",full.names=T)
#SampleConFile = read.table(RndConFiles[1])
#NumNodes      = tail(SampleConFile,n=1)[,2]
NumPairs      = (NumNodes-2)*(NumNodes-1)/2
AllCon.Vec    = GenerateAllConVec(NumNodes)


RndConWeightPIJList = GetPIJFromDir(RndConDir)
RndConWeight = RndConWeightPIJList$weight
RndConPIJ = RndConWeightPIJList$pij
RndConMat = RndConWeightPIJList$mat
NumFiles = length(RndConWeight)

NewConWeightPIJList = GetPIJFromDir(NewConDir_1)
NewConPIJ_1 = NewConWeightPIJList$pij

NewConWeightPIJList = GetPIJFromDir(NewConDir_2)
NewConPIJ_2 = NewConWeightPIJList$pij

#------------------------------------------------------------------
# Bootstrap for p-value
#------------------------------------------------------------------
BootstrapTimes = 1000
BootSampleSize = NumFiles
Bootstrap_1.Con = matrix(FALSE,nrow=NumPairs,ncol=BootstrapTimes)
Bootstrap_2.Con = matrix(FALSE,nrow=NumPairs,ncol=BootstrapTimes)

ptm <- proc.time()
for (i in 1:BootstrapTimes){
	if (i %% 100 == 0){
		cat(i,"/",BootstrapTimes,"\n"); flush.console();
	}
	Con.Samples.Ind = sample(1:NumFiles,BootSampleSize,replace=T)
	Con.pij = GetPIJFromMatWeight(RndConMat[,Con.Samples.Ind],RndConWeight[Con.Samples.Ind])
	Bootstrap_1.Con[,i] = (NewConPIJ_1 < Con.pij)
	Bootstrap_2.Con[,i] = (NewConPIJ_2 < Con.pij)
}
PVal_1.sim = apply(Bootstrap_1.Con,1,function(x){length(which(x))}) / BootstrapTimes
PVal_1.Vec = cbind(AllCon.Vec,PVal_1.sim)
PVal_2.sim = apply(Bootstrap_2.Con,1,function(x){length(which(x))}) / BootstrapTimes
PVal_2.Vec = cbind(AllCon.Vec,PVal_2.sim)
proc.time() - ptm

# ------------------------------------------------------------------------------------
# False Discovery Rate Resampling Adjustments 
# Yekutieli, D., & Benjamini, Y. (1999). Resampling-based false discovery rate controlling multiple test procedures for correlated test statistics. Journal of Statistical Planning and Inference, 82(1-2), 171-196. doi:10.1016/S0378-3758(99)00041-5
# Publishing, S. (2010). SAS/STAT 9. 22 User's Guide. The MIXED Procedure (Book Excerpt) (p. 228). SAS Institute.
# ------------------------------------------------------------------------------------

# different node space 
GetPValByFDR <- function(PVal.Vec){
	FDR_m = nrow(AllCon.Vec)
	Node.Space = unique(sort(AllCon.Vec[,2]-AllCon.Vec[,1]))

	New_P = rep(0,FDR_m)

	for (i in Node.Space){
		Node.Space.Ind = which(PVal.Vec[,2]-PVal.Vec[,1]==i,arr.ind=T)
		Node.Space.Size = length(Node.Space.Ind)
		if(Node.Space.Size==0){next}

		Node.Space.PVal = PVal.Vec[Node.Space.Ind,3]
		Node.Space.SortPVal = sort.int(PVal.Vec[Node.Space.Ind,3],index.return=T)
		Node.Space.SortPVal.ix = Node.Space.Ind[Node.Space.SortPVal$ix]
		Node.Space.SortPVal.x = Node.Space.SortPVal$x

		New_P_sub = rep(0,Node.Space.Size)
		New_P_sub[Node.Space.Size] = Node.Space.SortPVal.x[Node.Space.Size]
		if(Node.Space.Size>1){
			for(j in (Node.Space.Size-1):1){
				New_P_sub[j] = min(New_P_sub[j+1], (Node.Space.Size/j) * Node.Space.SortPVal.x[j] )
			}
		}

		New_P[Node.Space.SortPVal.ix] = New_P_sub
	}
	return(New_P)
}

PVal_1.FDR.Vec = cbind(AllCon.Vec,GetPValByFDR(PVal_1.Vec))
PVal_2.FDR.Vec = cbind(AllCon.Vec,GetPValByFDR(PVal_2.Vec))

Alpha = 0.01
AlphaStr = round(Alpha*100)
PickInd_1  = which(PVal_1.FDR.Vec[,3]<Alpha)
PickInd_2  = which(PVal_2.FDR.Vec[,3]<Alpha)

New.Vec_1 = AllCon.Vec[PickInd_1,]
New.Vec_2 = AllCon.Vec[PickInd_2,]
IdentifyVec <- function(Con.Vec,New.Vec){
	Match.Org.Vec = apply(Con.Vec,1,function(x){
			length(which(x[1] == New.Vec[,1] & x[2] == New.Vec[,2]))	
		}
	)
	Match.New.Vec = apply(New.Vec,1,function(x){
			length(which(x[1] == Con.Vec[,1] & x[2] == Con.Vec[,2]))	
		}
	)
	return(list(
		"common"=Con.Vec[Match.Org.Vec==1,],
		"new"=New.Vec[Match.New.Vec==0,],
		"missed"=Con.Vec[Match.Org.Vec==0,]))
}

VecList_1 = IdentifyVec(Con.Vec_1, New.Vec_1)
VecList_2 = IdentifyVec(Con.Vec_2, New.Vec_2)

WriteVecFile <- function(VecList,CELL){
	for (ItemName in names(VecList)){
		FN = paste(OutDir,"/",CELL,".",SampleSize,".c840.",ItemName,".a",AlphaStr,Type,".txt",sep="")
		cat(FN,"\n")
		write.table(file=FN,
			VecList[[ItemName]],quote=FALSE,row.names=FALSE,col.names=FALSE)
	}
}

OutDir = "tmp/"
WriteVecFile(VecList_1,"GM12878")
WriteVecFile(VecList_2,"K562")

