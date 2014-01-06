#!/usr/bin/R
###################################################################
# Get contact index according to segment (starting and end)
###################################################################

rm(list=ls())

library(reshape)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)

if (Sys.info()["sysname"] == "Darwin"){
  SampleSize    = "p1500_d3300"
}else{
  SampleSize    = args[1]
}

# 10 nm / 28 nm/kb = 357  bp \approx 350  bp
# BallSeqLen = 350
# 30 nm / 11 nm/kb =2727

PersistenceLength= as.numeric(sub("p(.*)_d.*","\\1",SampleSize))
if (PersistenceLength == 600){ # 10 nm diameter
  BallDiameter    = 100
  BallSeqLen      = 350
  SegmentBPLength = 2100
}else if (PersistenceLength == 900){ # 10 nm diameter
  BallDiameter    = 100
  BallSeqLen      = 350
  SegmentBPLength = 3200
}else if (PersistenceLength == 1200){ # 10 nm diameter
  BallDiameter    = 100
  BallSeqLen      = 350
  SegmentBPLength = 4300
}else if (PersistenceLength == 1700){ # 30 nm diameter
  BallDiameter    = 300
  BallSeqLen      = 2727
  SegmentBPLength = 15455
}else if (PersistenceLength == 1500){ # 30 nm diameter
  BallDiameter    = 300
  BallSeqLen      = 2727
  SegmentBPLength = 13636

}
#------------------------------------------------------------------
# Parameter
#------------------------------------------------------------------
SourcePath = "data/ENm008/"
# CELL = "GM12878"
ReferenceCountFile = paste("data/analysis/ENm008/",SampleSize,".comb.c.txt",sep="")
cat(ReferenceCountFile, "\n") 
DestDir = "tmp/"

#------------------------------------------------------------------
# Get MaxRead between GM12878 and K562
#------------------------------------------------------------------
ContactFrequencyFile   = paste(SourcePath,"Nature_2010_5CFrequencyCountsMatrix_ENm008_GM12878.txt",
    sep                = "")
ContactFrequencyData_1 = read.table(ContactFrequencyFile)

ContactFrequencyFile   = paste(SourcePath,"Nature_2010_5CFrequencyCountsMatrix_ENm008_K562.txt",
    sep                = "")
ContactFrequencyData_2 = read.table(ContactFrequencyFile)

#------------------------------------------------------------------
# Read Contact Frequency
#------------------------------------------------------------------
ContactFrequencyData = ContactFrequencyData_1

Row_Ind_Start        = as.numeric(sub(".*chr16:(.*)-(.*)","\\1",rownames(ContactFrequencyData)))
Row_Ind_End          = as.numeric(sub(".*chr16:(.*)-(.*)","\\2",rownames(ContactFrequencyData)))
Row_Ind_Diff         = Row_Ind_Start[2:(length(Row_Ind_Start))] - Row_Ind_End[1:(length(Row_Ind_End)-1)] 

Col_Ind_Start        = as.numeric(sub(".*chr16\\.(.*)\\.(.*)","\\1",colnames(ContactFrequencyData)))
Col_Ind_End          = as.numeric(sub(".*chr16\\.(.*)\\.(.*)","\\2",colnames(ContactFrequencyData)))
Col_Ind_Diff         = Col_Ind_Start[2:(length(Col_Ind_Start))] - Col_Ind_End[1:(length(Col_Ind_End)-1)] 



#------------------------------------------------------------------
# Reconstruct the reading count matrix, gap is considered.
#------------------------------------------------------------------
# Start from Col_Col_Start, Col_Ind_Start = 1
MinLength            = min(Row_Ind_End-Row_Ind_Start,Col_Ind_End - Col_Ind_Start)

Row_Vec              = cbind(Row_Ind_Start,Row_Ind_End)
Col_Vec              = cbind(Col_Ind_Start,Col_Ind_End)
All_Vec              = rbind(Row_Vec,Col_Vec)
All_Vec_Sort         = sort(All_Vec[,1],index.return=T)
All_Vec              = All_Vec[All_Vec_Sort$ix,]
All_Vec_Diff         = (All_Vec[,1])[2:nrow(All_Vec)] - (All_Vec[,2])[1:(nrow(All_Vec)-1)]

# expand the gap sequence length less than min seq length
MinInd               = which(All_Vec_Diff <MinLength & All_Vec_Diff != 0)
All_Vec[MinInd,2]    = All_Vec[MinInd,2] + All_Vec_Diff[MinInd]
All_Vec_Diff[MinInd] = 0

# insert the gap and resort, get segment interval
Gap_Vec              = cbind((All_Vec[,2])[1:(nrow(All_Vec)-1)], (All_Vec[,2])[1:(nrow(All_Vec)-1)] + All_Vec_Diff)
All_Vec              = rbind(All_Vec,Gap_Vec[which(All_Vec_Diff!=0),])
All_Vec_Sort         = sort(All_Vec[,1],index.return=T)
All_Vec              = All_Vec[All_Vec_Sort$ix,]

nSize                = nrow(All_Vec)
IntervalPoints       = unique(sort(c(All_Vec)))


#------------------------------
# Output Segemntindex NodeStart and NodeEnd
#------------------------------
# of_name = paste(DestDir,CELL,".SegInd.Node_StartEnd.txt",sep="")
# write.table(file=of_name,
#     cbind(1:nrow(All_Vec),All_Vec),quote=FALSE,row.names=FALSE,col.names=FALSE)

#----------------------------------------------------------------------------
# Map Experimental Index to equal segment index
#----------------------------------------------------------------------------
#------------------------------
# Persistence length
#------------------------------
# Dekker, J. (2008). Mapping in vivo chromatin interactions in yeast suggests an
# extended chromatin fiber with regional variation in compaction The Journal of
# biological chemistry, 283(50), 34532-34540. doi:10.1074/jbc.M806479200 the
# mass density 28 nm/kb the persistence length is 58-118 nm, the mean is 88 nm
#------------------------------
# Fussner, E., Ching, R. W., & Bazett-Jones, D. P. (2011). Living without 30nm
# chromatin fibers. Trends in biochemical sciences, 36(1), 1-6.
# doi:10.1016/j.tibs.2010.09.002 10 nm chromatin fiber
#------------------------------
# 90 nm / 28 nm/kb = 3214 bp \approx 3200 bp 10 nm / 28 nm/kb = 357  bp \approx
# 350  bp One segment length 88 nm - 10 nm = 78 nm, containing 3000-350 bp the
# diameter of ball is 10 nm. the ball contains 350 bp
#------------------------------

EndPoint = max(IntervalPoints)
NumSegments = ceiling(EndPoint/SegmentBPLength)
SegmentsInterval = c(1,c(1:NumSegments)*SegmentBPLength)

SegStartEnd.Forward.List = list()
for (i in 1:nrow(All_Vec)){
    SegStartEndInd = findInterval(All_Vec[i,], SegmentsInterval)
    SegStartEnd.Forward.List[[i]] <- c(SegStartEndInd[1] : SegStartEndInd[2])
}

SegStartEnd.Reverse.List = list()
for (i in 1:NumSegments) {
    SegStartEnd.Reverse.List[[i]] = seq_along(SegStartEnd.Forward.List)[sapply(SegStartEnd.Forward.List, FUN=function(X) i %in% X)] 
}




#----------------------------------------------------------------------------
# Get new contact table
# pos1 pos2 count
#----------------------------------------------------------------------------
MapRawData <- function(CFData){
  # remove small length segment
  # only some row index, because
  IgnoredRow_Ind    = which(Row_Ind_End-Row_Ind_Start <= BallSeqLen)
  NewRow_Ind        = setdiff(1:length(Row_Ind_End),IgnoredRow_Ind)
  # column index do not change
  IgnoredCol_Ind    = which(Col_Ind_End-Col_Ind_Start <= BallSeqLen)
  NewCol_Ind        = setdiff(1:length(Col_Ind_Start),IgnoredCol_Ind)

  REVERSE_POS       = Col_Ind_End[NewCol_Ind]
  REVERSE_POS_END   = Col_Ind_End[NewCol_Ind]
  REVERSE_POS_START = Col_Ind_Start[NewCol_Ind]

  FORWARD_POS       = Row_Ind_End[NewRow_Ind]
  FORWARD_POS_START = Row_Ind_Start[NewRow_Ind]
  FORWARD_POS_END   = Row_Ind_End[NewRow_Ind]

  CFData            = CFData[NewRow_Ind,NewCol_Ind]

  SparseData = NULL
  for (i in 1:length(FORWARD_POS)){
    for (j in 1:length(REVERSE_POS)){
      ipos        = FORWARD_POS[i]
      ipos_start  = FORWARD_POS_START[i]
      ipos_end    = FORWARD_POS_END[i]
      jpos        = REVERSE_POS[j]
      jpos_start  = REVERSE_POS_START[j]
      jpos_end    = REVERSE_POS_END[j]
      # remove neighbor node
      if (ipos_end == jpos_start || ipos_start == jpos_end) {
        cat(
          paste("F",which(Row_Ind_End == FORWARD_POS_END[i]),"-","R",which(Col_Ind_End == REVERSE_POS_END[j]),sep=""),"\t",
          paste(FORWARD_POS_END[i],"-",REVERSE_POS_END[j],sep=""),"\t",CFData[i,j],"\n")
        next;
      }
      # if (ipos_start == jpos_end) {cat("F",i,"-",FORWARD_POS_END[i],"R",j,"-",REVERSE_POS_END[j],"\t",CFData[i,j],"\n"); next;}
      if (ipos < jpos){
        SparseData = rbind(SparseData,c(ipos,jpos,CFData[i,j]))
      }else{
        SparseData = rbind(SparseData,c(jpos,ipos,CFData[i,j]))

      }
    }
  }

  # sort first and second
  SparseData = SparseData[sort.int(SparseData[,1],index.return=T)$ix,]
  SparseData = SparseData[sort.int(SparseData[,2],index.return=T)$ix,]
  
  return (SparseData)
}



MapRawData_PureRaw <- function(CFData){
  # remove small length segment
  # only some row index, because
  # IgnoredRow_Ind    = which(Row_Ind_End-Row_Ind_Start <= BallSeqLen)
  # NewRow_Ind        = setdiff(1:length(Row_Ind_End),IgnoredRow_Ind)
  # # column index do not change
  # IgnoredCol_Ind    = which(Col_Ind_End-Col_Ind_Start <= BallSeqLen)
  # NewCol_Ind        = setdiff(1:length(Col_Ind_Start),IgnoredCol_Ind)

	NewRow_Ind = 1:length(Row_Ind_End)
	NewCol_Ind = 1:length(Col_Ind_Start)

  REVERSE_POS       = Col_Ind_End[NewCol_Ind]
  REVERSE_POS_END   = Col_Ind_End[NewCol_Ind]
  REVERSE_POS_START = Col_Ind_Start[NewCol_Ind]

  FORWARD_POS       = Row_Ind_End[NewRow_Ind]
  FORWARD_POS_START = Row_Ind_Start[NewRow_Ind]
  FORWARD_POS_END   = Row_Ind_End[NewRow_Ind]

  CFData            = CFData[NewRow_Ind,NewCol_Ind]

  SparseData = NULL
  for (i in 1:length(FORWARD_POS)){
    for (j in 1:length(REVERSE_POS)){
      ipos        = FORWARD_POS[i]
      ipos_start  = FORWARD_POS_START[i]
      ipos_end    = FORWARD_POS_END[i]
      jpos        = REVERSE_POS[j]
      jpos_start  = REVERSE_POS_START[j]
      jpos_end    = REVERSE_POS_END[j]
      # remove neighbor node
      # if (ipos_end == jpos_start || ipos_start == jpos_end) {
      #   cat(
      #     paste("F",which(Row_Ind_End == FORWARD_POS_END[i]),"-","R",which(Col_Ind_End == REVERSE_POS_END[j]),sep=""),"\t",
      #     paste(FORWARD_POS_END[i],"-",REVERSE_POS_END[j],sep=""),"\t",CFData[i,j],"\n")
      #   next;
      # }
      # if (ipos_start == jpos_end) {cat("F",i,"-",FORWARD_POS_END[i],"R",j,"-",REVERSE_POS_END[j],"\t",CFData[i,j],"\n"); next;}
      if (ipos < jpos){
        SparseData = rbind(SparseData,c(ipos,jpos,CFData[i,j]))
      }else{
        SparseData = rbind(SparseData,c(jpos,ipos,CFData[i,j]))

      }
    }
  }

  # sort first and second
  SparseData = SparseData[sort.int(SparseData[,1],index.return=T)$ix,]
  SparseData = SparseData[sort.int(SparseData[,2],index.return=T)$ix,]
  
  return (SparseData)
}

# CFData_1 = MapRawData_PureRaw(ContactFrequencyData_1)
# CFData_2 = MapRawData_PureRaw(ContactFrequencyData_2)
WriteCircos_RawConnection <- function(FileName,Data){

	LinkStr = ""
	ii = 1
	for (iCount in 1:nrow(Data)){
		if (Data[iCount,3]!=0) {
			LinkStr = paste (LinkStr, sprintf("segdup%04d hs16 ", ii), Data[iCount,1], " ", Data[iCount,1], " id=2\n",sep="")
			LinkStr = paste (LinkStr, sprintf("segdup%04d hs16 ", ii), Data[iCount,2], " ", Data[iCount,2], " id=2\n",sep="")
			ii = ii + 1
		}
	}
	
	write.table(LinkStr, file=FileName,quote=F,row.names=F,col.names=F)
}
# 
# FN_1 = "/tmp/GM12878.Raw.Connection.Circos.txt"
# WriteCircos_RawConnection(FN_1, CFData_1)
# FN_2 = "/tmp/K562.Raw.Connection.Circos.txt"
# WriteCircos_RawConnection(FN_2, CFData_2)



CFData_1 = MapRawData(ContactFrequencyData_1)
CFData_2 = MapRawData(ContactFrequencyData_2)

#----------------------------------------------------------------------------
# Get Index cut by persistence length
Points_Pos = unique(sort(CFData_1[,1:2]))
Points_Pos = c(1,Points_Pos)
Len_Diff = Points_Pos[2:length(Points_Pos)] - Points_Pos[1:(length(Points_Pos)-1)] 

Seg_bp = NULL
for (i in 1:length(Len_Diff)){
  Len_Cur = Len_Diff[i]
  while(Len_Cur >0){
    if (Len_Cur < SegmentBPLength){
      if (Len_Cur < BallSeqLen){
        Seg_bp[length(Seg_bp)] = Seg_bp[length(Seg_bp)] + Len_Cur
      }else{
        Seg_bp = c(Seg_bp,Len_Cur)
      }
    }else{
      Seg_bp = c(Seg_bp, SegmentBPLength)
    }
    Len_Cur = Len_Cur - SegmentBPLength
  }
}
Seg_Len = Seg_bp * (PersistenceLength/SegmentBPLength)
All_Points_Pos = cumsum(c(1,Seg_bp))
Points_Map_Ind = which(cumsum(c(1,Seg_bp)) %in% Points_Pos)

# of_name = paste(DestDir,"ENm008.","p",PersistenceLength,".equ.seg.len.txt",sep="")
# write.table(file=of_name,sprintf("%.3f",Seg_Len),row.names=FALSE,col.names=FALSE,quote=FALSE)


# transform exp point index to equavalent segment point index
GetEquPointIndFromExpPointInd <- function(ExpPointInd){
  return(Points_Map_Ind[ExpPointInd])
}

GetEquPointIndFromExpPointPos <- function(ExpPointPos){
  ExpPointInd = NULL
  for (i in 1:length(ExpPointPos)){
    ExpPointInd = c(ExpPointInd,which(Points_Pos == ExpPointPos[i]))
  }
  return (GetEquPointIndFromExpPointInd(ExpPointInd))
}

GetExpPointIndFromEquPointInd <- function(EquPointInd){
  ExpPointInd = NULL
  for (i in 1:length(EquPointInd)){
    ExpPointInd = c(ExpPointInd,which(Points_Map_Ind==EquPointInd[i]))
  }
  return(ExpPointInd)
}

GetExpPointPosFromEquPointInd <- function(EquPointInd){
  ExpPointPos = NULL
  for (i in 1:length(EquPointInd)){
    ExpPointPos = c(ExpPointPos,
      Points_Pos[GetExpPointIndFromEquPointInd(EquPointInd[i])])
  }
  return(ExpPointPos)
}

GetClosestEquPointIndFromExpPointPos <- function(ExpPointPos){
	ExpPointInd = NULL
	for (i in 1:length(ExpPointPos)){
		Points_Diff = abs(All_Points_Pos - ExpPointPos[i])
		ExpPointInd = c(ExpPointInd,which(Points_Diff == min(Points_Diff)))
	}
	return (ExpPointInd)
}

Equ_Contact_Point_Ind = cbind(
  GetEquPointIndFromExpPointPos(CFData_1[,1]),
  GetEquPointIndFromExpPointPos(CFData_1[,2]))
  
# remove equ neighbor segment
EquSegNonNeighborInd = which(Equ_Contact_Point_Ind[,2]-Equ_Contact_Point_Ind[,1]!=1)
Equ_Contact_Point_Ind = Equ_Contact_Point_Ind[EquSegNonNeighborInd,]

CFData_1 = CFData_1[EquSegNonNeighborInd,]
CFData_2 = CFData_2[EquSegNonNeighborInd,]

CF_1.Vec = cbind(Equ_Contact_Point_Ind,CFData_1[,3])
CF_2.Vec = cbind(Equ_Contact_Point_Ind,CFData_2[,3])
###################################################################
#  Node set
NumNodes = length(All_Points_Pos)
# IgnoreNodes = c(1,11,16,19,20,28,30,31,39,44,45,51)
# NodeSet = setdiff(1:NumNodes,IgnoreNodes)
NodeSet = 1:NumNodes

###################################################################
# read dist file to get neighbour numbers
DataPath = "/Users/yunxu/workspace/projects/chromatin/result/analysis/ENm008/chain/weight.dist/par.1_1_1_1/"

NodeNodeList = t(combn(1:NumNodes,2))
NodeNodeList = NodeNodeList[which(NodeNodeList[,2]-NodeNodeList[,1]!=1),]

GetContactNumMedianSd <- function(Cell){
	ContactNum.sts = matrix(0,nrow=NumNodes,ncol=2)
	for (iNode in 1:NumNodes){
		cat(iNode,"\n");
		flush.console()
		iNodeInd = which(NodeNodeList[,1]==iNode | NodeNodeList[,2]==iNode)
		NodeNodeFileList = paste(DataPath,"/",Cell,".p1500_d3300/",apply(NodeNodeList[iNodeInd,],1,function(x){paste(x,collapse="_")}),".txt",sep="")

		SampleLen = nrow(read.table(NodeNodeFileList[1]))
		SampleWeight = read.table(NodeNodeFileList[1])[,2]

		ContactMat = matrix(0, nrow=SampleLen, ncol=length(NodeNodeFileList))
		for (j in 1:length(NodeNodeFileList)){
			jFile = NodeNodeFileList[j]
			DistTmp = read.table(jFile)[,1]

			ContactMat[which(DistTmp < 840),j] = 1
		}
		ContactNum = rep(apply(ContactMat,1,sum) ,SampleWeight)
		ContactNum.sts[iNode,] = c(median(ContactNum),sd(ContactNum))
	}
	return(ContactNum.sts)
}

# GM_ContactNum.sts = GetContactNumMedianSd("GM12878")
# K_ContactNum.sts = GetContactNumMedianSd("K562")


###################################################################
#  read Contact Index file come from significant
sts_path = "../result/analysis/ENm008/chain/sts/pij/"
ContactIndexFile_GM = paste(sts_path,"GM12878.p1500_d3300.c840.con.pij.all.contactindex.txt",sep="")
ContactIndexFile_K = paste(sts_path,"K562.p1500_d3300.c840.con.pij.all.contactindex.txt",sep="")
#PijFile_GM = paste(sts_path,"GM12878.p1500_d3300.c840.con.pij.triplet.allall.txt",sep="")
#PijFile_K = paste(sts_path,"K562.p1500_d3300.c840.con.pij.triplet.allall.txt",sep="")

# read pij come from significant contact
PijFile_GM_Significant =
  paste(sts_path,"GM12878.p1500_d3300.c840.con.pij.all.triplet.txt",sep="")
PijFile_K_Significant =
  paste(sts_path,"K562.p1500_d3300.c840.con.pij.all.triplet.txt",sep="")

CI_GM = read.table(ContactIndexFile_GM)[,2]
CI_K = read.table(ContactIndexFile_K)[,2]

Pij_GM_Significant = read.table(PijFile_GM_Significant)
Pij_K_Significant  = read.table(PijFile_K_Significant)

Pij_GM = Pij_GM_Significant
Pij_K = Pij_K_Significant
###################################################################
# read referennce file
Con.Vec            = CF_1.Vec[,1:2]
SampleSize    = "p1500_d3300"
ReferenceCountFile = paste("data/analysis/ENm008/",SampleSize,".comb.c.txt",sep="")
RF.Vec = read.table(ReferenceCountFile)
# check if same contact in same row
for ( i in 1:nrow(CF_1.Vec)){
  MatchInd = which(CF_1.Vec[i,1]==RF.Vec[,1] & CF_1.Vec[i,2]==RF.Vec[,2])
  if(length(MatchInd) == 0 ){
    stop(paste("contact do not match in ",i," row", sep=" "))
  }
}
# replace zero with minimal
Contact.RF.Count.Vec = RF.Vec[,3:ncol(RF.Vec)]

Contact.RF.Count.Vec =  apply(Contact.RF.Count.Vec,2,function(x)
      { matchind = which(x==0)
        if(length(matchind)){x[matchind]= min(x[which(x!=0)])}
        return(x) })

RF.Vec = cbind(Con.Vec, Contact.RF.Count.Vec)
RF.Vec.Scale = cbind(Con.Vec, apply(Contact.RF.Count.Vec,2,function(x){x/sum(x)}))


MaxNode = 54
GetSumNodeProb <- function(ScaledVec){
  # ScaledVec = RF.Vec.Scale
  SumNodeProb = NULL
  for (Nodei in 1:MaxNode){
    SumNodeProb = c(SumNodeProb,sum(ScaledVec[which(ScaledVec[,1]==Nodei | ScaledVec[,2] ==Nodei),3]))
  }
  return (SumNodeProb)
}

RandSumNodeProb = GetSumNodeProb(RF.Vec.Scale)
###################################################################
#  read Dnase | CTCF

circos_path = "data/ENm008/circos/"

Read_MarkerEnrichment <- function(file){
	# file = CTCFFile_K
	Data_Raw = read.table(file)[,2:4]
	MiddlePoints = apply(Data_Raw,1,function(x){mean(x[1:2])})
	MiddlePointsMapToNodeIndex = GetClosestEquPointIndFromExpPointPos(MiddlePoints)

	Marker_Enrichment = NULL
	for (i in 1:NumNodes){
		MiddlePoints_Ind = which(MiddlePointsMapToNodeIndex == i)
		if(length(MiddlePoints_Ind) != 0){
			Marker_Enrichment = rbind(Marker_Enrichment, c(i,max(Data_Raw[MiddlePoints_Ind,3])))
		}
		else{
			Marker_Enrichment = rbind(Marker_Enrichment, c(i,0))
		}
	}
	return (Marker_Enrichment)
}

# ----------------------------------------------------------------------------------------
DNaseFile_Duke_GM     =  paste(circos_path , "wgEncodeDukeDNaseSeqPeaksGm12878V3_circos.txt"        , sep="")
DNaseFile_Duke_K      =  paste(circos_path , "wgEncodeDukeDNaseSeqPeaksK562V3_circos.txt"           , sep="")
DNaseFile_UW_GM       =  paste(circos_path , "wgEncodeUwDnaseSeqPeaksRep1Gm12878_circos.txt"        , sep="")
DNaseFile_UW_K        =  paste(circos_path , "wgEncodeUwDnaseSeqPeaksRep1K562_circos.txt"           , sep="")

CTCFFile_Broad_GM     =  paste(circos_path , "wgEncodeBroadChipSeqPeaksGm12878Ctcf_circos.txt"      , sep="")
CTCFFile_Broad_K      =  paste(circos_path , "wgEncodeBroadChipSeqPeaksK562Ctcf_circos.txt"         , sep="")
CTCFFile_UW_GM        =  paste(circos_path , "wgEncodeUwChIPSeqPeaksRep1Gm12878Ctcf_circos.txt"     , sep="")
CTCFFile_UW_K         =  paste(circos_path , "wgEncodeUwChIPSeqPeaksRep1K562Ctcf_circos.txt"        , sep="")
CTCFFile_Uta_GM       =  paste(circos_path , "wgEncodeUtaChIPseqPeaksGm12878CtcfV3_circos.txt"      , sep="")
CTCFFile_Uta_K        =  paste(circos_path , "wgEncodeUtaChIPseqPeaksK562CtcfV3_circos.txt"         , sep="")

H3K4me1File_Broad_GM  =  paste(circos_path , "wgEncodeBroadChipSeqPeaksGm12878H3k4me1_circos.txt"   , sep="")
H3K4me1File_Broad_K   =  paste(circos_path , "wgEncodeBroadChipSeqPeaksK562H3k4me1_circos.txt"      , sep="")

H3K4me2File_Broad_GM  =  paste(circos_path , "wgEncodeBroadChipSeqPeaksGm12878H3k4me2_circos.txt"   , sep="")
H3K4me2File_Broad_K   =  paste(circos_path , "wgEncodeBroadChipSeqPeaksK562H3k4me2_circos.txt"      , sep="")

H3K4me3File_Broad_GM  =  paste(circos_path , "wgEncodeBroadChipSeqPeaksGm12878H3k4me3_circos.txt"   , sep="")
H3K4me3File_Broad_K   =  paste(circos_path , "wgEncodeBroadChipSeqPeaksK562H3k4me3_circos.txt"      , sep="")
H3K4me3File_UW_GM     =  paste(circos_path , "wgEncodeUwChIPSeqPeaksRep1Gm12878H3k4me3_circos.txt"  , sep="")
H3K4me3File_UW_K      =  paste(circos_path , "wgEncodeUwChIPSeqPeaksRep1K562H3k4me3_circos.txt"     , sep="")

H3K9acFile_Broad_GM   =  paste(circos_path , "wgEncodeBroadChipSeqPeaksGm12878H3k9ac_circos.txt"    , sep="")
H3K9acFile_Broad_K    =  paste(circos_path , "wgEncodeBroadChipSeqPeaksK562H3k9ac_circos.txt"       , sep="")

H3K27me3File_Broad_GM =  paste(circos_path , "wgEncodeBroadChipSeqPeaksGm12878H3k27me3_circos.txt"  , sep="")
H3K27me3File_Broad_K  =  paste(circos_path , "wgEncodeBroadChipSeqPeaksK562H3k27me3_circos.txt"     , sep="")
H3K27me3File_UW_GM    =  paste(circos_path , "wgEncodeUwChIPSeqPeaksRep1Gm12878H3k27me3_circos.txt" , sep="")
H3K27me3File_UW_K     =  paste(circos_path , "wgEncodeUwChIPSeqPeaksRep1K562H3k27me3_circos.txt"    , sep="")

H3K27acFile_Broad_GM  =  paste(circos_path , "wgEncodeBroadChipSeqPeaksGm12878H3k27ac_circos.txt"   , sep="")
H3K27acFile_Broad_K   =  paste(circos_path , "wgEncodeBroadChipSeqPeaksK562H3k27ac_circos.txt"      , sep="")

H3K36me3File_Broad_GM =  paste(circos_path , "wgEncodeBroadChipSeqPeaksGm12878H3k36me3_circos.txt"  , sep="")
H3K36me3File_Broad_K  =  paste(circos_path , "wgEncodeBroadChipSeqPeaksK562H3k36me3_circos.txt"     , sep="")
H3K36me3File_UW_GM    =  paste(circos_path , "wgEncodeUwChIPSeqPeaksRep1Gm12878H3k36me3_circos.txt" , sep="")
H3K36me3File_UW_K     =  paste(circos_path , "wgEncodeUwChIPSeqPeaksRep1K562H3k36me3_circos.txt"    , sep="")

H4K20me1File_Broad_GM =  paste(circos_path , "wgEncodeBroadChipSeqPeaksGm12878H4k20me1_circos.txt"  , sep="")
H4K20me1File_Broad_K  =  paste(circos_path , "wgEncodeBroadChipSeqPeaksK562H4k20me1_circos.txt"     , sep="")

# ----------------------------------------------------------------------------------------
DNase_Duke_GM     = Read_MarkerEnrichment ( DNaseFile_Duke_GM     ) [,2]
DNase_Duke_K      = Read_MarkerEnrichment ( DNaseFile_Duke_K      ) [,2]
DNase_UW_GM       = Read_MarkerEnrichment ( DNaseFile_UW_GM       ) [,2]
DNase_UW_K        = Read_MarkerEnrichment ( DNaseFile_UW_K        ) [,2]

CTCF_Broad_GM     = Read_MarkerEnrichment ( CTCFFile_Broad_GM     ) [,2]
CTCF_Broad_K      = Read_MarkerEnrichment ( CTCFFile_Broad_K      ) [,2]
CTCF_UW_GM        = Read_MarkerEnrichment ( CTCFFile_UW_GM        ) [,2]
CTCF_UW_K         = Read_MarkerEnrichment ( CTCFFile_UW_K         ) [,2]
CTCF_Uta_GM       = Read_MarkerEnrichment ( CTCFFile_Uta_GM       ) [,2]
CTCF_Uta_K        = Read_MarkerEnrichment ( CTCFFile_Uta_K        ) [,2]

H3K4me1_Broad_GM  = Read_MarkerEnrichment ( H3K4me1File_Broad_GM  ) [,2]
H3K4me1_Broad_K   = Read_MarkerEnrichment ( H3K4me1File_Broad_K   ) [,2]

H3K4me2_Broad_GM  = Read_MarkerEnrichment ( H3K4me2File_Broad_GM  ) [,2]
H3K4me2_Broad_K   = Read_MarkerEnrichment ( H3K4me2File_Broad_K   ) [,2]

H3K4me3_Broad_GM  = Read_MarkerEnrichment ( H3K4me3File_Broad_GM  ) [,2]
H3K4me3_Broad_K   = Read_MarkerEnrichment ( H3K4me3File_Broad_K   ) [,2]
H3K4me3_UW_GM     = Read_MarkerEnrichment ( H3K4me3File_UW_GM     ) [,2]
H3K4me3_UW_K      = Read_MarkerEnrichment ( H3K4me3File_UW_K      ) [,2]

H3K9ac_Broad_GM   = Read_MarkerEnrichment ( H3K9acFile_Broad_GM   ) [,2]
H3K9ac_Broad_K    = Read_MarkerEnrichment ( H3K9acFile_Broad_K    ) [,2]

H3K27me3_Broad_GM = Read_MarkerEnrichment ( H3K27me3File_Broad_GM ) [,2]
H3K27me3_Broad_K  = Read_MarkerEnrichment ( H3K27me3File_Broad_K  ) [,2]
H3K27me3_UW_GM    = Read_MarkerEnrichment ( H3K27me3File_UW_GM    ) [,2]
H3K27me3_UW_K     = Read_MarkerEnrichment ( H3K27me3File_UW_K     ) [,2]

H3K27ac_Broad_GM  = Read_MarkerEnrichment ( H3K27acFile_Broad_GM  ) [,2]
H3K27ac_Broad_K   = Read_MarkerEnrichment ( H3K27acFile_Broad_K   ) [,2]

H3K36me3_Broad_GM = Read_MarkerEnrichment ( H3K36me3File_Broad_GM ) [,2]
H3K36me3_Broad_K  = Read_MarkerEnrichment ( H3K36me3File_Broad_K  ) [,2]
H3K36me3_UW_GM    = Read_MarkerEnrichment ( H3K36me3File_UW_GM    ) [,2]
H3K36me3_UW_K     = Read_MarkerEnrichment ( H3K36me3File_UW_K     ) [,2]

H4K20me1_Broad_GM = Read_MarkerEnrichment ( H4K20me1File_Broad_GM ) [,2]
H4K20me1_Broad_K  = Read_MarkerEnrichment ( H4K20me1File_Broad_K  ) [,2]

# ----------------------------------------------------------------------------------------
DNase_Duke_GM_Normalize     = DNase_Duke_GM     / max(DNase_Duke_GM)
DNase_Duke_K_Normalize      = DNase_Duke_K      / max(DNase_Duke_K)
DNase_UW_GM_Normalize       = DNase_UW_GM       / max(DNase_UW_GM)
DNase_UW_K_Normalize        = DNase_UW_K        / max(DNase_UW_K)

CTCF_Broad_GM_Normalize     = CTCF_Broad_GM     / max(CTCF_Broad_GM)
CTCF_Broad_K_Normalize      = CTCF_Broad_K      / max(CTCF_Broad_K)
CTCF_UW_GM_Normalize        = CTCF_UW_GM        / max(CTCF_UW_GM)
CTCF_UW_K_Normalize         = CTCF_UW_K         / max(CTCF_UW_K)
CTCF_Uta_GM_Normalize       = CTCF_Uta_GM       / max(CTCF_Uta_GM)
CTCF_Uta_K_Normalize        = CTCF_Uta_K        / max(CTCF_Uta_K)

H3K4me1_Broad_GM_Normalize  = H3K4me1_Broad_GM  / max(H3K4me1_Broad_GM)
H3K4me1_Broad_K_Normalize   = H3K4me1_Broad_K   / max(H3K4me1_Broad_K)

H3K4me2_Broad_GM_Normalize  = H3K4me2_Broad_GM  / max(H3K4me2_Broad_GM)
H3K4me2_Broad_K_Normalize   = H3K4me2_Broad_K   / max(H3K4me2_Broad_K)

H3K4me3_Broad_GM_Normalize  = H3K4me3_Broad_GM  / max(H3K4me3_Broad_GM)
H3K4me3_Broad_K_Normalize   = H3K4me3_Broad_K   / max(H3K4me3_Broad_K)
H3K4me3_UW_GM_Normalize     = H3K4me3_UW_GM     / max(H3K4me3_UW_GM)
H3K4me3_UW_K_Normalize      = H3K4me3_UW_K      / max(H3K4me3_UW_K)

H3K9ac_Broad_GM_Normalize   = H3K9ac_Broad_GM   / max(H3K9ac_Broad_GM)
H3K9ac_Broad_K_Normalize    = H3K9ac_Broad_K    / max(H3K9ac_Broad_K)

H3K27me3_Broad_GM_Normalize = H3K27me3_Broad_GM / max(H3K27me3_Broad_GM)
H3K27me3_Broad_K_Normalize  = H3K27me3_Broad_K  / max(H3K27me3_Broad_K)
H3K27me3_UW_GM_Normalize    = H3K27me3_UW_GM    / max(H3K27me3_UW_GM)
H3K27me3_UW_K_Normalize     = H3K27me3_UW_K     / max(H3K27me3_UW_K)

H3K27ac_Broad_GM_Normalize  = H3K27ac_Broad_GM  / max(H3K27ac_Broad_GM)
H3K27ac_Broad_K_Normalize   = H3K27ac_Broad_K   / max(H3K27ac_Broad_K)

H3K36me3_Broad_GM_Normalize = H3K36me3_Broad_GM / max(H3K36me3_Broad_GM)
H3K36me3_Broad_K_Normalize  = H3K36me3_Broad_K  / max(H3K36me3_Broad_K)
H3K36me3_UW_GM_Normalize    = H3K36me3_UW_GM    / max(H3K36me3_UW_GM)
H3K36me3_UW_K_Normalize     = H3K36me3_UW_K     / max(H3K36me3_UW_K)

H4K20me1_Broad_GM_Normalize = H4K20me1_Broad_GM / max(H4K20me1_Broad_GM)
H4K20me1_Broad_K_Normalize  = H4K20me1_Broad_K  / max(H4K20me1_Broad_K)

# ----------------------------------------------------------------------------------------
# dev.new ( width=12,height=4                    )
# plot    ( DNase_Duke_GM_Normalize,type="o"     )
# lines   ( DNase_UW_GM_Normalize,col="red"      )
# 
# dev.new ( width=12,height=4                    )
# plot    ( DNase_Duke_K_Normalize,type="o"      )
# lines   ( DNase_UW_K_Normalize,col="red"       )
# 
# dev.new ( width=12,height=4                    )
# plot    ( CTCF_Broad_GM_Normalize,type="o"     )
# lines   ( CTCF_UW_GM_Normalize,col="red"       )
# lines   ( CTCF_Uta_GM_Normalize,col="blue"     )
# 
# dev.new ( width=12,height=4                    )
# plot    ( CTCF_Broad_K_Normalize,type="o"      )
# lines   ( CTCF_UW_K_Normalize,col="red"        )
# lines   ( CTCF_Uta_K_Normalize,col="blue"      )
# 
# dev.new ( width=12,height=4                    )
# plot    ( H3K4me3_Broad_GM_Normalize,type="o"  )
# lines   ( H3K4me3_UW_GM_Normalize,col="red"    )
# 
# dev.new ( width=12,height=4                    )
# plot    ( H3K4me3_Broad_K_Normalize,type="o"   )
# lines   ( H3K4me3_UW_K_Normalize,col="red"     )
# 
# dev.new ( width=12,height=4                    )
# plot    ( H3K27me3_Broad_GM_Normalize,type="o" )
# lines   ( H3K27me3_UW_GM_Normalize,col="red"   )
# 
# dev.new ( width=12,height=4                    )
# plot    ( H3K27me3_Broad_K_Normalize,type="o"  )
# lines   ( H3K27me3_UW_K_Normalize,col="red"    )
# 
# dev.new ( width=12,height=4                    )
# plot    ( H3K36me3_Broad_GM_Normalize,type="o" )
# lines   ( H3K36me3_UW_GM_Normalize,col="red"   )
# 
# dev.new ( width=12,height=4                    )
# plot    ( H3K36me3_Broad_K_Normalize,type="o"  )
# lines   ( H3K36me3_UW_K_Normalize,col="red"    )
# 
# dev.new ( width=12,height=4                    )
# plot    ( H3K4me1_Broad_GM_Normalize,type="o"  )
# lines   ( H3K4me2_Broad_GM_Normalize,col="red" )
#
# dev.new ( width=12,height=4                    )
# plot    ( H3K4me1_Broad_K_Normalize,type="o"  )
# lines   ( H3K4me2_Broad_K_Normalize,col="red" )
#----------------------------------------------------------------------------
GetGeneNameAndRange <- function(){
  # All_Points_Pos = c(1,5693,15091,18344,29756,44231,50868,55911,64056,74448,88084,95256,100530,104687,109838,123474,131220,134334,147970,161606,167103,171769,185994,189074,203353,217802,225341,238977,247214,260850,274486,277942,289198,303867,310528,314538,327240,334889,348525,352385,360924,372670,380160,393796,407432,418222,421291,433293,445126,454365,468001,483412,496513,499411)
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

# ----------------------------------------------------------------------------------------
# ggplot2
# one figure for GM and K comparison
ggplot_Enrichment_oneplot_GMandK <- function(df_long, ylab="") {

  MyPal = brewer.pal(9,"Set1")
  dodge = position_dodge(width=0.5)

  gp <- ggplot(df_long, aes(x= Node, y=value, color= variable,group=variable))
  gp = gp +
    theme_bw() +
    theme(
      text                 = element_text(size  = 15),
      axis.text.x          = element_text(angle = 30, hjust = 1,vjust = 1),
      legend.justification = c(1,1),
      legend.position      = c(1,.9) ) +
    labs(
      x= paste("Primer sites"),
      y = ylab) +
    scale_x_continuous(breaks=1:NumNodes,labels=xlabel) +
    scale_color_manual("",
     values = c( "GM12878" = MyPal[1], "K562" = MyPal[2]))  
  gp = gp + geom_linerange(
    subset=.(variable %in% c("GM12878","K562")),
    aes(ymin=0,ymax=value,color=variable),
    position= dodge,size=2)

}

ggplot_draw <- function(gp){
  dev.new(width=16,height=4)
  gp
}

# ----------------------------------------------------------------------------------------
# CTCF analysis                                                                          # 
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# CTCF raw GM and K data
df = data.frame(
    Node    = 1:NumNodes,
    GM12878 = CTCF_Broad_GM,
    K562    = CTCF_Broad_K)
df_long = melt(df,id="Node")
gp = ggplot_Enrichment_oneplot_GMandK(df_long,ylab="CTCF enrichment")
ggplot_draw(gp)
# ----------------------------------------------------------------------------------------
# H3K4me1 raw GM and K data
df = data.frame(
    Node    = 1:NumNodes,
    GM12878 = H3K4me1_Broad_GM,
    K562    = H3K4me1_Broad_K)
df_long = melt(df,id="Node")
gp = ggplot_Enrichment_oneplot_GMandK(df_long,ylab="H3K4me1 enrichment")
ggplot_draw(gp)
# ----------------------------------------------------------------------------------------
# CTCF normalize GM and K data
df = data.frame(
    Node    = 1:NumNodes,
    GM12878 = CTCF_Broad_GM_Normalize,
    K562    = CTCF_Broad_K_Normalize)
df_long = melt(df,id="Node")
gp = ggplot_Enrichment_oneplot_GMandK(df_long,ylab="CTCF normalized enrichment")
ggplot_draw(gp)


# ----------------------------------------------------------------------------------------
GetNumCTCFInteractions <- function(Cell, Node_Ind){
	# Cell = "GM12878"
	# Node_Ind = CTCF_GM_nozeroind

	Node_Ind = sort(Node_Ind)
	NNList = t(combn(Node_Ind,2))
	NNList = NNList[NNList[,2]-NNList[,1]!=1,]

	NNFileList = paste(DataPath,"/",Cell,".p1500_d3300/",apply(NNList,1,function(x){paste(x,collapse="_")}),".txt",sep="")
	SampleLen = nrow(read.table(NNFileList[1]))
	SampleWeight = read.table(NNFileList[1])[,2]
	
	ContactMat = matrix(0,nrow=SampleLen,ncol=nrow(NNList))
	for (j in 1:length(NNFileList)){
		jFile = NNFileList[j]
		DistTmp = read.table(jFile)[,1]
		ContactMat[which(DistTmp< 840),j] = 1
		cat(j, "/", length(NNFileList), "\n")
		flush.console()
	}
	return(mean(rep(apply(ContactMat,1,sum) ,SampleWeight)))
}

GetProbCTCFInteractions <- function(NumExp, Node_Ind){
	# NumExp = NumCTCFInteraction_GM
	# Node_Ind = CTCF_GM_nozeroind

	Node_Ind = sort(Node_Ind)
	NNList = t(combn(Node_Ind,2))
	NNList = NNList[NNList[,2]-NNList[,1]!=1,]
	DataPath_Rand = "/Users/yunxu/workspace/projects/chromatin/result/analysis/ENm008/ensemble/weight.dist"

	NNFileList = paste(DataPath_Rand,"/",apply(NNList,1,function(x){paste(x,collapse="_")}),".txt",sep="")
	SampleLen = nrow(read.table(NNFileList[1]))
	SampleWeight = read.table(NNFileList[1])[,2]
	
	ContactMat = matrix(0,nrow=SampleLen,ncol=nrow(NNList))
	for (j in 1:length(NNFileList)){
		jFile = NNFileList[j]
		DistTmp = read.table(jFile)[,1]
		ContactMat[which(DistTmp< 860),j] = 1
		cat(j, "/", length(NNFileList), "\n")
		flush.console()
	}

	NumBoots = 100
	MeanNum = NULL
	for (i in 1:NumBoots){
		SampleInd = sample(1:SampleLen,SampleLen,replace=T)
		MeanNum = c(MeanNum, mean(rep(apply(ContactMat[SampleInd,],1,sum) ,SampleWeight[SampleInd])) )
		cat(i, "/", NumBoots, "\n")
		flush.console()
	}
	
	Prob = length(which(MeanNum > NumExp))/NumBoots
	
	return(Prob)
}

# ----------------------------------------------------------------------------------------
CTCF_Broad_GM_nozeroind = which(CTCF_Broad_GM_Normalize != 0)
CTCF_Broad_K_nozeroind = which(CTCF_Broad_K_Normalize != 0)
# NumCTCFInteraction_Broad_GM = GetNumCTCFInteractions("GM12878",CTCF_Broad_GM_nozeroind)
# NumCTCFInteraction_Broad_K = GetNumCTCFInteractions("K562",CTCF_Broad_K_nozeroind)
# cat("The average number of CTCF node interaction for GM is \t",
#     NumCTCFInteraction_Broad_GM)
# cat("The average number of CTCF node interaction for K is \t",
#         NumCTCFInteraction_Broad_K)

# Prob_CTCF_Broad_GM =
# GetProbCTCFInteractions(NumCTCFInteraction_Broad_GM,CTCF_Broad_GM_nozeroind)
# Prob_CTCF_Broad_K =
# GetProbCTCFInteractions(NumCTCFInteraction_Broad_K,CTCF_Broad_K_nozeroind)

# ----------------------------------------------------------------------------------------
CTCF_Pij_Broad_GM = Pij_GM[
  which(Pij_GM[,1] %in% CTCF_Broad_GM_nozeroind & 
        Pij_GM[,2] %in% CTCF_Broad_GM_nozeroind),]
CTCF_Pij_Broad_K = Pij_K[
  which(Pij_K[,1] %in% CTCF_Broad_K_nozeroind & 
        Pij_K[,2] %in% CTCF_Broad_K_nozeroind),]

# ----------------------------------------------------------------------------------------
# GetCTCFCorrelation <- function(Pij_prd, CTCF_Exp, iNodeVec=NULL){
# #  Pij_prd  = Pij_GM
# #  CTCF_Exp = CTCF_Broad_GM
# #  iNodeVec = NULL
# 
#   CTCF_threshold = 20
#   CTCF_threshold_ind = which (CTCF_Exp > CTCF_threshold)
#   
#   # remove node 1 for its inconsistence
#   if (CTCF_threshold_ind[1] == 1){
#     CTCF_threshold_ind = CTCF_threshold_ind[-1]  
#   }
#   
#   if (length(iNodeVec) != 0){
#     iNodeList = iNodeVec
#   } else {
#     iNodeList = CTCF_threshold_ind
#   }
# 
#   CTCF_List = list()
# #  print(substitute(Pij_prd))
#   for (i in 1:length(iNodeList)){
#     iNode = iNodeList[i]
# 
#     Pij_prd_iNode = Pij_prd[Pij_prd[,1] == iNode | Pij_prd[,2] == iNode,]
#     Pij_prd_iNode_CTCF =
#      Pij_prd_iNode[Pij_prd_iNode[,1] %in% CTCF_threshold_ind &
#              Pij_prd_iNode[,2] %in% CTCF_threshold_ind,]
# 
# 
#     if(nrow(Pij_prd_iNode_CTCF)!=0){
#       CTCF_enr = apply(Pij_prd_iNode_CTCF,1,
#             function(x){
#               CTCF_Exp[x[1]] * CTCF_Exp[x[2]]
#             })
#       cat(iNode, cor(CTCF_enr, Pij_prd_iNode_CTCF[,3]),"\n")
# 
#     } else {
#       cat (iNode, "has no significant contact\n")
#       CTCF_enr = NULL
#     }
#     CTCF_List [[as.character(iNode)]] = cbind(Pij_prd_iNode_CTCF, CTCF_enr)
#     print (CTCF_List [[as.character(iNode)]])
#     cat("--------------\n")
#   }
#   return(CTCF_List)
# }

#CTCF_Broad_GM_Pij_Exp = GetCTCFCorrelation(Pij_GM, CTCF_Broad_GM)
#CTCF_Broad_K_Pij_Exp = GetCTCFCorrelation(Pij_K, CTCF_Broad_K)
# -----------------------------------------------------------------------------------
# calculate correlation between CTCF enrichment and Pij
# GetCTCFCorrelation_ByPij<- function(Pij_prd, CTCF_Exp, CTCF_threshold, Pij_threshold){
# #  Pij_prd  = Pij_K
# #  CTCF_Exp = CTCF_Broad_K
# #  CTCF_threshold = 15
# #  Pij_threshold = .4
# 
#   CTCF_threshold_ind = which (CTCF_Exp > CTCF_threshold)
#   if (CTCF_threshold_ind[1] == 1){
#     CTCF_threshold_ind = CTCF_threshold_ind[-1]
#   }
# #  print(CTCF_threshold_ind)
#   Pij_prd_threshold = Pij_prd[Pij_prd[,3]>Pij_threshold,]
#   Pij_prd_CTCF = 
#     Pij_prd_threshold[ 
#       Pij_prd_threshold[,1] %in% CTCF_threshold_ind &
#       Pij_prd_threshold[,2] %in% CTCF_threshold_ind ,]
#   
#   CTCF_enr = apply(Pij_prd_CTCF,1,
#             function(x){
#               CTCF_Exp[x[1]] * CTCF_Exp[x[2]]
#             })
#   cat(
#     CTCF_threshold, 
#     Pij_threshold,
#     sprintf("%2.3f",cor(Pij_prd_CTCF[,3],CTCF_enr)),
#     length(CTCF_threshold_ind),
#     length(CTCF_enr),
# #    nrow(Pij_prd_threshold)
#     "\n",sep="  "
#     )
# #  print(Pij_prd_CTCF)
# #  return(cor(Pij_prd_CTCF[,3],CTCF_enr))
#   return(Pij_prd_CTCF)
# 
# }
# 
# CTCF_threshold_list = seq(0,20,by=5)
# Pij_threshold_list = seq(0.3,0.5,by=.1)
# # CTCF_threshold_list = seq(0,40,by=5)
# # Pij_threshold_list = seq(0.2,1.0,by=.1)
# #CTCF_threshold_list = 15
# #Pij_threshold_list = 0.4
# 
# for (CTCF_threshold in CTCF_threshold_list){
#   for (Pij_threshold in Pij_threshold_list){
#     GetCTCFCorrelation_ByPij(Pij_GM, CTCF_Broad_GM, CTCF_threshold, Pij_threshold)
#   }
# }
# 
# for (CTCF_threshold in CTCF_threshold_list){
#   for (Pij_threshold in Pij_threshold_list){
#     GetCTCFCorrelation_ByPij(Pij_K, CTCF_Broad_K, CTCF_threshold, Pij_threshold)
#   }
# }
# ----------------------------------------------------------------------------------------
# GetCTCFCorrelation_ByRnd<- function(Pij_prd, CTCF_Exp, CTCF_threshold, Pij_threshold){
# 
# #  Pij_prd  = Pij_GM
# #  CTCF_Exp = CTCF_Broad_GM
# #  Pij_prd  = Pij_K
# #  CTCF_Exp = CTCF_Broad_K
# #  CTCF_threshold = 15
# #  Pij_threshold = .4
# 
#   CTCF_threshold_ind = which (CTCF_Exp > CTCF_threshold)
#   Pij_prd_threshold = Pij_prd[Pij_prd[,3]>Pij_threshold,]
#   Pij_prd_CTCF = 
#     Pij_prd_threshold[ 
#       Pij_prd_threshold[,1] %in% CTCF_threshold_ind &
#       Pij_prd_threshold[,2] %in% CTCF_threshold_ind ,]
# 
#   CTCF_Ind = unique(sort(unlist(Pij_prd_CTCF[,1:2])))
# 
#   CTCF_enr = apply(Pij_prd_CTCF,1,
#             function(x){
#               CTCF_Exp[x[1]] * CTCF_Exp[x[2]]
#             })
#   Cor_Real = cor(Pij_prd_CTCF[,3],CTCF_enr) 
# 
#   Cor_List = NULL
#   BootTimes = 100000
# 
#   Sign_Greater_Vec = NULL
# #  Node_Selected = setdiff(1:NumNodes, CTCF_threshold_ind)
#   Node_Selected = 1:NumNodes
#   for (i in 1:BootTimes){
#     Sign_Greater = 0
#     if(i %% (BootTimes/10) == 0){
#       cat(i,"/", BootTimes,"\n")
#       flush.console()
#     }
#     Pij_Rnd = as.data.frame(NULL)
#     Rnd_Ind = sample(Node_Selected,length(CTCF_threshold_ind ))
#     Enr_Rnd = sample(CTCF_Exp[CTCF_threshold_ind],length(CTCF_threshold_ind ))
#     Pij_Rnd = Pij_prd_threshold[
#         Pij_prd_threshold[,1] %in% Rnd_Ind &
#         Pij_prd_threshold[,2] %in% Rnd_Ind ,]
#     if (nrow(Pij_Rnd) >= nrow(Pij_prd_CTCF)){
#       CTCF_Rnd = apply(Pij_Rnd,1,function(x){
#         Enr_Rnd[which(Rnd_Ind==x[1])]* Enr_Rnd[which(Rnd_Ind==x[2])]
#         })
#       Cor_Rnd = cor(CTCF_Rnd, Pij_Rnd[,3])
#       if (Cor_Rnd > Cor_Real){
#         Sign_Greater = 1      
#       }
#     }
#     Sign_Greater_Vec = c(Sign_Greater_Vec, Sign_Greater)    
#   }
#   P_Value = sum(Sign_Greater_Vec)/BootTimes
# 
#   return(P_Value)
# }


# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------
# GetProbCTCF_ByRndEnsemble<- function(Pij_prd, CTCF_Exp, CTCF_threshold, Pij_threshold){
# 
# #  Pij_prd  = Pij_GM
# #  CTCF_Exp = CTCF_Broad_GM
# #  CTCF_threshold = 15
# #  Pij_threshold = .4
# 
#   CTCF_threshold_ind = which (CTCF_Exp > CTCF_threshold)
#   if (CTCF_threshold_ind[1] == 1){
#     CTCF_threshold_ind = CTCF_threshold_ind[-1]
#   }
#   print(CTCF_threshold_ind)
# 
#   Pij_prd_threshold = Pij_prd[Pij_prd[,3]>Pij_threshold,]
#   Pij_prd_CTCF = 
#     Pij_prd_threshold[ 
#       Pij_prd_threshold[,1] %in% CTCF_threshold_ind &
#       Pij_prd_threshold[,2] %in% CTCF_threshold_ind ,]
#   
#   CTCF_enr = apply(Pij_prd_CTCF,1,
#             function(x){
#               CTCF_Exp[x[1]] * CTCF_Exp[x[2]]
#             })
#   Cor_Real = cor(Pij_prd_CTCF[,3] , CTCF_enr)
# 
#   Pair_Vec = t(combn(CTCF_threshold_ind,2))
#   Pair_Vec = Pair_Vec[Pair_Vec[,2] - Pair_Vec[,1] != 1,]
# 
# 	DataPath_Rand = "/Users/yunxu/workspace/projects/chromatin/result/analysis/ENm008/ensemble/weight.dist"
#   PairFileList = paste(DataPath_Rand,"/",apply(Pair_Vec,1,function(x){paste(x,collapse="_")}),".txt",sep="")
# 
#   SampleFile = read.table(PairFileList[1])
#   NRow = nrow(SampleFile)
#   Weight_Vec = SampleFile[,2]
# 
#   DistMat = matrix(0,nrow=NRow, ncol=nrow(Pair_Vec))
#   for (i in 1:length(PairFileList)){
#    PairFile = PairFileList[i]
#    DistMat[,i] = read.table(PairFile)[,1]
#   }
# 
#   SignMat = matrix(0,nrow=NRow, ncol=nrow(Pair_Vec))
#   SignMat[which(DistMat < 840)] = 1
# 
#   CTCF_Pair = apply(Pair_Vec,1,function(x){
#                    CTCF_Exp[x[1]] * CTCF_Exp[x[2]]
#             })
#  
#   BootTimes = 10000
#   Num_Vec = rep(0, BootTimes)
#   NumSample = round(nrow(SignMat)/10)
#   for (i in 1:BootTimes){
#     if (i %% (BootTimes/10) == 0) {
#       cat (i ,"/", BootTimes, "\n")
#       flush.console()
#     }
# 
#     SampleInd = sample(1:nrow(SignMat),NumSample)
#     Pij_Sample = apply(SignMat[SampleInd,] * Weight_Vec[SampleInd], 2, sum) / sum(Weight_Vec[SampleInd])
#     
#     Pij_greater_Ind = which(Pij_Sample > Pij_threshold)
#     if (length(Pij_greater_Ind) >= nrow(Pij_prd_CTCF)){
#       cat(length(Pij_greater_Ind),"\n")
#       if (cor(Pij_Sample[Pij_greater_Ind], CTCF_Pair[Pij_greater_Ind]) > Cor_Real){
#         Num_Vec[i] = 1
#       } 
#     }
# 
#   }
# 
#   P_Value = sum(Num_Vec)/BootTimes
# 
#   return(P_Value)
# 
# }
# # ----------------------------------------------------------------------------------------
# CTCF_threshold = 15
# Pij_threshold = 0.4
# 
# (CTCF_Cor_GM = GetCTCFCorrelation_ByPij(Pij_GM, CTCF_Broad_GM, CTCF_threshold, Pij_threshold))
# (CTCF_Cor_GM_Pval = GetProbCTCF_ByRndEnsemble(Pij_GM, CTCF_Broad_GM, CTCF_threshold, Pij_threshold))
# 
# (CTCF_Cor_K = GetCTCFCorrelation_ByPij(Pij_K, CTCF_Broad_K, CTCF_threshold, Pij_threshold))
# (CTCF_Cor_K_Pval = GetProbCTCF_ByRndEnsemble(Pij_K, CTCF_Broad_K, CTCF_threshold, Pij_threshold))






# ----------------------------------------------------------------------------------------
GetStatCTCF_ByRndEnsemble<- function(Pij_prd, CTCF_Exp, CTCF_threshold, Pij_threshold){
#  Pij_prd  = Pij_GM
#  CTCF_Exp = CTCF_Broad_GM
#  Pij_prd  = Pij_K
#  CTCF_Exp = CTCF_Broad_K
#  CTCF_threshold = 0
#  Pij_threshold = .3

  CTCF_threshold_ind = which (CTCF_Exp > CTCF_threshold)
  if (CTCF_threshold_ind[1] == 1){
    CTCF_threshold_ind = CTCF_threshold_ind[-1]
  }
  print(CTCF_threshold_ind)

  Pij_prd_threshold = Pij_prd[Pij_prd[,3]>Pij_threshold,]
  Pij_prd_CTCF = 
    Pij_prd_threshold[ 
      Pij_prd_threshold[,1] %in% CTCF_threshold_ind &
      Pij_prd_threshold[,2] %in% CTCF_threshold_ind ,]
  

  Pair_Vec = t(combn(CTCF_threshold_ind,2))
  Pair_Vec = Pair_Vec[Pair_Vec[,2] - Pair_Vec[,1] != 1,]

	DataPath_Rand = "/Users/yunxu/workspace/projects/chromatin/result/analysis/ENm008/ensemble/weight.dist"
  PairFileList = paste(DataPath_Rand,"/",apply(Pair_Vec,1,function(x){paste(x,collapse="_")}),".txt",sep="")

  SampleFile = read.table(PairFileList[1])
  NRow = nrow(SampleFile)
  Weight_Vec = SampleFile[,2]

  DistMat = matrix(0,nrow=NRow, ncol=nrow(Pair_Vec))
  for (i in 1:length(PairFileList)){
   PairFile = PairFileList[i]
   DistMat[,i] = read.table(PairFile)[,1]
  }

  SignMat = matrix(0,nrow=NRow, ncol=nrow(Pair_Vec))
  SignMat[which(DistMat < 840)] = 1

  Num_CTCFCTCF = rep(apply(SignMat,1, sum), Weight_Vec)
  PValue = length(Num_CTCFCTCF[Num_CTCFCTCF > nrow(Pij_prd_CTCF)])/ length(Num_CTCFCTCF)
  cat(CTCF_threshold, Pij_threshold, nrow(Pij_prd_CTCF), PValue, "\n", sep=" ")
  return(list(PijCTCFPair=Pij_prd_CTCF,NumCTCFPair = nrow(Pij_prd_CTCF),NumCTCFPairVec=Num_CTCFCTCF, PValue=PValue))
}

# ----------------------------------------------------------------------------------------
# Get pvalue from random chain and draw it
CTCF_threshold_list = seq(0,20,by=5)
Pij_threshold_list = seq(0.3,0.5,by=.1)
CTCF_threshold = 0
Pij_threshold = 0.3

cat ("GM-------------------------\n")
StatCTCF_GM = GetStatCTCF_ByRndEnsemble(Pij_GM, CTCF_Broad_GM, CTCF_threshold, Pij_threshold)
StatCTCF_GM$PijCTCFPair
cat ("K-------------------------\n")
StatCTCF_K = GetStatCTCF_ByRndEnsemble(Pij_K, CTCF_Broad_K, CTCF_threshold, Pij_threshold)
StatCTCF_K$PijCTCFPair

Draw_CTCFPair_Hist <- function(StatCTCF_Cell){

  MyPal = brewer.pal(9,"Set1")
  CellName = unlist(strsplit(as.character(substitute(StatCTCF_Cell)),"_"))[2]
  if (CellName=="GM"){
    CellColor = MyPal[1]
  } else{
    CellColor = MyPal[2]
  }

  df = data.frame(Cell=CellName,NumCTCFPair=StatCTCF_Cell$NumCTCFPairVec)
  
  gp = ggplot(df,aes(x=NumCTCFPair)) 
  gp = gp +
    theme_bw() +
    theme(text = element_text(size  = 15)) +
    labs( x= paste("Number of CTCF pairs"), y = "Count") +
    scale_x_continuous(limits=c(0,30)) + 
    scale_y_continuous(limits=c(0,260000))
  gp = gp + geom_histogram(colour = "darkgreen", fill = "white",binwidth=0.5)
  gp = gp + geom_vline(xintercept = StatCTCF_Cell$NumCTCFPair, color=CellColor,size = 1)
  gp = gp + annotate("text", 
         label = paste( "P-Value =", sprintf("%.2e",StatCTCF_Cell$PValue)),
         x     = StatCTCF_Cell$NumCTCFPair,
         y     = max(table(StatCTCF_Cell$NumCTCFPairVec)) *5/6)
  return(gp)
}

gp_hist_GM = Draw_CTCFPair_Hist(StatCTCF_GM)
gp_hist_K = Draw_CTCFPair_Hist(StatCTCF_K)
grid.arrange(gp_hist_GM,gp_hist_K)

# ----------------------------------------------------------------------------------------
# output for circos
Get_CTCF_CircosStr <- function(StatCTCF_Cell){
#  StatCTCF_Cell = StatCTCF_GM
  PijCTCFPair = StatCTCF_Cell$PijCTCFPair
  row.names(PijCTCFPair) = 1:nrow(PijCTCFPair)

  idx = 1
  str_circos = apply(PijCTCFPair,1,function(x){
    str_head = sprintf("segdup%04d hs16",idx);
    str1 = paste(str_head,
      All_Points_Pos[x[1]],All_Points_Pos[x[1]], 
      sprintf("pij=%.2f",x[3]))
    str2 = paste(str_head,
      All_Points_Pos[x[2]],All_Points_Pos[x[2]], 
      sprintf("pij=%.2f",x[3]))
         
    idx <<- idx + 1;
    paste(str1,str2,sep="\n")
    }
    )
  write.table(paste(str_circos,collapse="\n"),col.names=F,row.names=F,quote=F)
}
# GM12878.CTCFPair.Pij.Connection.txt
# K562.CTCFPair.Pij.Connection.txt
Get_CTCF_CircosStr(StatCTCF_GM)
Get_CTCF_CircosStr(StatCTCF_K)

# ----------------------------------------------------------------------------------------
# get CTCF one pair p-value
Get_CTCF_PairPValue <- function(){
  DataPlain = "8  11  0.2747
  8 22  0.1999
  8 40  NA
  9 11  0.3319
  9 22  NA
  11  12  NA
  11  22  0.2002
  11  40  NA
  12  22  0.5052
  13  22  0.5291
  15  22  0.3417
  19  20  NA
  43  52  0.3615
  "
  PijNodePair = read.table(textConnection(DataPlain))
  PijNodePair = PijNodePair[!is.na(PijNodePair[,3]) & PijNodePair[,3]>0.3,]
  print(PijNodePair)

	DataPath_Rand = "/Users/yunxu/workspace/projects/chromatin/result/analysis/ENm008/ensemble/weight.dist"
  PairFileList = paste(DataPath_Rand,"/",apply(PijNodePair,1,function(x){paste(x[1:2],collapse="_")}),".txt",sep="")

  SampleFile = read.table(PairFileList[1])
  NRow = nrow(SampleFile)
  Weight_Vec = SampleFile[,2]

  DistMat = matrix(0,nrow=NRow, ncol=nrow(PijNodePair))
  for (i in 1:length(PairFileList)){
   PairFile = PairFileList[i]
   DistMat[,i] = read.table(PairFile)[,1]
  }

  SignMat = matrix(0,nrow=NRow, ncol=nrow(PijNodePair))
  SignMat[which(DistMat < 840)] = 1

  BootTimes = 1000
  NumSample = NRow

  Res_Mat = matrix(0,nrow=BootTimes,ncol=nrow(PijNodePair))
  for (i in 1:BootTimes){
    if (i %% (BootTimes/10) == 0){
      cat (i , "/", BootTimes, "\n")
      flush.console()
    }
    Rnd_Ind = sample(1:NRow,NumSample,replace=T)
    Rnd_Sample = SignMat[Rnd_Ind,] * Weight_Vec[Rnd_Ind] 
    Res_Mat[i,] = apply(Rnd_Sample,2,sum)/sum(Weight_Vec)
  }

  idx = 1
  apply(Res_Mat,2,function(x){
    PValue = length(which(x > PijNodePair[idx,3])) / BootTimes
    idx <<- idx + 1
    PValue
    }) 

}
# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------
# output for circos
Get_CHIAPET_CTCF_CircosStr <- function(){
  DataPlain = "8  11  0.2747
  8 22  0.1999
  8 40  NA
  9 11  0.3319
  9 22  NA
  11  12  NA
  11  22  0.2002
  11  40  NA
  12  22  0.5052
  13  22  0.5291
  15  22  0.3417
  19  20  NA
  43  52  0.3615
  "
  PijNodePair = read.table(textConnection(DataPlain))

  idx = 1
  str_circos = apply(PijNodePair,1,function(x){
    if(is.na(x[3])){
      x[3] = -1
    } 
    if(x[2]-x[1] == 1){
      x[3] = 1
    }
    str_head = sprintf("segdup%04d hs16",idx);
    str1 = paste(str_head,
      All_Points_Pos[x[1]],All_Points_Pos[x[1]], 
      sprintf("pij=%.2f",x[3]))
    str2 = paste(str_head,
      All_Points_Pos[x[2]],All_Points_Pos[x[2]], 
      sprintf("pij=%.2f",x[3]))
         
    idx <<- idx + 1;
    paste(str1,str2,sep="\n")
    }
    )
  write.table(paste(str_circos,collapse="\n"),col.names=F,row.names=F,quote=F)
}
# ----------------------------------------------------------------------------------------
# quit R
q()















#  Pij_prd  = Pij_GM
#  CTCF_Exp = CTCF_Broad_GM
#  Pij_prd  = Pij_K
#  CTCF_Exp = CTCF_Broad_K
#  CTCF_threshold = 0
#  Pij_threshold = .4


# GetNumCTCFBinding_ByPij<- function(Pij_prd, CTCF_Exp, CTCF_threshold_arg=0, Pij_threshold_arg=0.3){
# 
# #  Pij_prd  = Pij_GM
# #  CTCF_Exp = CTCF_Broad_GM
# #  CTCF_threshold_arg = 0
# #  Pij_threshold_arg = .3
# 
#   CTCF_threshold_ind = which (CTCF_Exp > CTCF_threshold_arg)
# 
#   Pij_prd_threshold = Pij_prd[Pij_prd[,3]>Pij_threshold_arg,]
#   Pij_prd_CTCF = 
#     Pij_prd_threshold[ 
#       Pij_prd_threshold[,1] %in% CTCF_threshold_ind &
#       Pij_prd_threshold[,2] %in% CTCF_threshold_ind ,]
#   return(nrow(Pij_prd_CTCF))  
# }
# 
# GetNumCTCFBinding_ByRand<- function(Pij_prd, CTCF_Exp, CTCF_threshold_arg=0, Pij_threshold_arg=0.3){
# 
#   Pij_prd  = Pij_K
#   CTCF_Exp = CTCF_Broad_K
#   CTCF_threshold_arg = 0
#   Pij_threshold_arg = .3
# 
#   CTCF_threshold_ind = which (CTCF_Exp > CTCF_threshold_arg)
#   Pij_prd_threshold = Pij_prd[Pij_prd[,3]>Pij_threshold_arg,]
# 
#   NumCTCF = length(CTCF_threshold_ind)
#   BootTimes = 10000
# 
#   NumList = NULL
#   for (i in 1:BootTimes){
#     if (i %% 100 ==0){
#       cat (i,"/",BootTimes,"\n")
#       flush.console()
#     }
#     Rnd_Ind = sort(sample(1:NumNodes,NumCTCF))
# 
#     Pij_prd_CTCF_Rnd = 
#       Pij_prd_threshold[ 
#         Pij_prd_threshold[,1] %in% Rnd_Ind &
#         Pij_prd_threshold[,2] %in% Rnd_Ind ,]
#     NumList = c(NumList, nrow(Pij_prd_CTCF_Rnd))
#   }
#   return(NumList)  
# }
# aa = GetNumCTCFBinding_ByPij(Pij_GM, CTCF_Broad_GM, 0, 0.3)

# ----------------------------------------------------------------------------------------
Draw_CTCF_Pij_Exp <- function(CTCF_Exp, CTCF_Cor_Cell, Pij_Cell){

#  CTCF_Exp = CTCF_Broad_GM
#  CTCF_Cor_Cell = CTCF_Cor_GM
#  Pij_Cell = Pij_GM

  CTCF_threshold_ind = which (CTCF_Exp > CTCF_threshold)
  if (CTCF_threshold_ind[1] == 1){
        CTCF_threshold_ind = CTCF_threshold_ind[-1]
  }

  iNodeList = unique(sort(CTCF_Cor_Cell[,1]))

  par(mfrow=c(length(iNodeList),1),mar=c(2, 4, 0, 1))
  for (iNode in iNodeList){
    Pij_Subset = Pij_Cell[ which(Pij_Cell[,1] == iNode | Pij_Cell[,2] == iNode), ]
    Pij_x = apply(Pij_Subset,1,function(x){if(x[1]==iNode){x[2]}else{x[1]} })
    Pij_y = Pij_Subset[,3]

    CTCF_ind = which(
          (Pij_Subset[,2] %in% CTCF_threshold_ind & Pij_Subset[,1] == iNode) |
          (Pij_Subset[,1] %in% CTCF_threshold_ind & Pij_Subset[,2] == iNode)  )
    CTCF_x = Pij_Subset[CTCF_ind,1:2]
    CTCF_x = unique(sort(unlist(CTCF_x)))
    CTCF_x = CTCF_x[(CTCF_x)!=iNode]
    CTCF_y = Pij_Subset[CTCF_ind,3]

    Cor_vec = t(apply(CTCF_Cor_Cell[which(CTCF_Cor_Cell[,1]==iNode | CTCF_Cor_Cell[,2]==iNode),,drop=F],1,
      function(x){
        if(x[1] == iNode){ 
          return (x[c(2,3)])
        } 
        if(x[2] == iNode){
          return (x[c(1,3)])
        }
          }))
#    dev.new(width=16,height=4)
    plot(Pij_x,Pij_y,type="o",ylim=c(0,1), xlim=c(1,NumNodes))
    points(CTCF_x, CTCF_y, col="red", pch=16)
    points(Cor_vec, col="blue", pch=5,cex=2)
    abline(v=iNode,col="pink",lwd=4)
    xTick = unique(sort(c(seq(0,54,by=10),CTCF_threshold_ind)))
    axis(side=1,at=xTick,labels=xTick)
  }

}
  
dev.new(width=7.7,height= 8.7*5/10)
Draw_CTCF_Pij_Exp(CTCF_Broad_GM,CTCF_Cor_GM,Pij_GM)
dev.new(width=7.7,height= 8.7)
Draw_CTCF_Pij_Exp(CTCF_Broad_K,CTCF_Cor_K,Pij_K)

# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------
# For H3K4me1
GetSign <- function(vec){
	sign_vec = NULL
	for (i in 1:length(vec)){
		if (i == 1){
			sign = 1
		}else{
			if ( vec[i] >= vec[i-1]  ){
				sign = 1
			}else{
				sign = -1
			}
		}
		sign_vec = c(sign_vec, sign)
	}
	return(sign_vec)
}

CheckSameTrend <- function(diff_ind, diff_sign){
	# diff_ind = DNase_Duke_GM_nozeroind
	# diff_sign = DNase_Duke_GM_sign - CI_GM_sign
	zero_ind = which(diff_sign == 0)
	TwoEndVec = matrix(1,nrow=1,ncol=2)
	for (i in 2:length(zero_ind)){
		if (zero_ind[i]-zero_ind[i-1] == 1){
			TwoEndVec[nrow(TwoEndVec),2] = zero_ind[i]
		} else {
			TwoEndVec = rbind(TwoEndVec,c(zero_ind[i],zero_ind[i]))
		}
	}
	ContinousNum = 3
	TwoEndVec = TwoEndVec[which(TwoEndVec[,2]-TwoEndVec[,1]+1 >= ContinousNum),,drop=FALSE]

	if(length(TwoEndVec)==0){
		NodesEndVec = NULL
	} else {
		NodesEndVec = t(apply(TwoEndVec,1,
			function(x){
				if(x[1]==1){return(c(diff_ind[ x[1] ], diff_ind[ x[2] ]))
			}else{
				return(c(diff_ind[ x[1]-1 ], diff_ind[ x[2] ]))}
			}))
		
	}
	return(NodesEndVec)
}

GetRegionCorrelation <- function(Enrichment,ContactIndex,SameTrend){
	# Enrichment = DNase_UW_GM_Normalize
	# ContactIndex = CI_GM_Normalize
	# SameTrend = DNase_UW_GM_sametrend

	CorVec = rep(0,nrow(SameTrend))
	NonZeroInd = which(Enrichment != 0 )
	for (i in 1:nrow(SameTrend)){
		RegionInd = NonZeroInd[which(NonZeroInd >= SameTrend[i,1] & NonZeroInd <= SameTrend[i,2])]
		CorVec[i] = try(cor(Enrichment[RegionInd],ContactIndex[RegionInd]))
		if(is.na(CorVec[i])){
			CorVec[i] = "NA"
		}
	}
	return(CorVec)
}
# ----------------------------------------------------------------------------------------
DNase_Duke_GM_nozeroind = which(DNase_Duke_GM_Normalize != 0)
DNase_Duke_GM_sign = GetSign(DNase_Duke_GM_Normalize[DNase_Duke_GM_nozeroind])
CI_GM_sign = GetSign(CI_GM_Normalize[DNase_Duke_GM_nozeroind])
rbind(DNase_Duke_GM_nozeroind, DNase_Duke_GM_sign - CI_GM_sign)
DNase_Duke_GM_sametrend = CheckSameTrend(DNase_Duke_GM_nozeroind,DNase_Duke_GM_sign - CI_GM_sign)
DNase_Duke_GM_Cor = GetRegionCorrelation(DNase_Duke_GM_Normalize,CI_GM_Normalize,DNase_Duke_GM_sametrend)
# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------
# GetCTCFCorrelation_Max <- function(Pij_prd, CTCF_Exp){
# #  Pij_prd  = Pij_K
# #  CTCF_Exp = CTCF_Broad_K
#   CTCF_threshold = 20
#   CTCF_threshold_ind = which (CTCF_Exp > CTCF_threshold)
#   
#   # remove node 1 for its inconsistence
#   if (CTCF_threshold_ind[1] == 1){
#     CTCF_threshold_ind = CTCF_threshold_ind[-1]  
#   }
#   print(substitute(Pij_prd))
# 
#   Pij_hit_Max = NULL
#   for (i in 1:length(CTCF_threshold_ind)){
#     iNode = CTCF_threshold_ind [i]
# 
#     Pij_hit = NULL
#     CTCF_hit = NULL
#     for (j in 1:length(CTCF_threshold_ind)){
#       jNode = CTCF_threshold_ind[j]
#       # remove neighbour and self
#       if (iNode == jNode | abs(iNode-jNode) == 1){
#         next
#       }
#       
#       if (iNode < jNode){
#         Pij_hit = rbind(Pij_hit,c(jNode, 
#            Pij_prd[which(Pij_prd[,1] == iNode & Pij_prd[,2] == jNode),3] ))
#       } else {
#         Pij_hit = rbind(Pij_hit,c(jNode, 
#            Pij_prd[which(Pij_prd[,1] == jNode & Pij_prd[,2] == iNode),3] ))
#       }
#       CTCF_hit = rbind(CTCF_hit,c(jNode, CTCF_Exp[iNode] * CTCF_Exp[jNode]))
#     }
# #    cat ("\n")
# #    cat ( iNode, cor(CTCF_hit[,2],Pij_hit[,2]),"\n")
#     
#     Pij_CTCF_List = Pij_hit[which(Pij_hit[,1] > iNode),2]
#     if (length(Pij_CTCF_List)){
#       Pij_hit_Max = rbind(Pij_hit_Max,
#             c(iNode, Pij_hit[which(Pij_hit[,2] ==
#                              max(Pij_CTCF_List)),]))
#     }
#    # print(Pij_hit)
#   }
#   print (Pij_hit_Max)
#   CTCF_hit_Max =
#   cbind(Pij_hit_Max[,1:2],apply(Pij_hit_Max,1,function(x){CTCF_Exp[x[1]]*CTCF_Exp[x[2]]}))
#   print (CTCF_hit_Max)
#   print (cor(CTCF_hit_Max[,3], Pij_hit_Max[,3]))
# }
# GetCTCFCorrelation_Max(Pij_GM, CTCF_Broad_GM)
# GetCTCFCorrelation_Max(Pij_K, CTCF_Broad_K)
# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------






























# ----------------------------------------------------------------------------------------

###################################################################

GetPval <- function(NodeSet_GM, Marker_GM, CI_GM, NodeSet_K, Marker_K, CI_K){
	# NodeSet_GM = GM_Retain.Ind
	# NodeSet_K  = K_Retain.Ind
	# Marker_GM = DNase_Duke_GM
	# Marker_K  = DNase_Duke_K
	# Marker_GM = H3K4me1_GM
	# Marker_K  = H3K4me1_K
	
	Marker_GM.Sel = Marker_GM[NodeSet_GM]
	CI_GM.Sel = CI_GM[NodeSet_GM]
	Marker_K.Sel = Marker_K[NodeSet_K]
	CI_K.Sel = CI_K[NodeSet_K]
	
	NumPermute = 1000
	Correlation_GM = cor(Marker_GM.Sel, CI_GM.Sel)
	Correlation_K = cor(Marker_K.Sel, CI_K.Sel)

	CorMat = matrix(0, nrow=NumPermute, ncol=2)
	for (i in 1:NumPermute){
		NewInd_GM = sample(1:length(NodeSet_GM),length(NodeSet_GM))
		NewCor_GM = cor(Marker_GM.Sel[NewInd_GM], CI_GM.Sel)
		NewInd_K = sample(1:length(NodeSet_K),length(NodeSet_K))
		NewCor_K = cor(Marker_K.Sel[NewInd_K], CI_K.Sel)
		CorMat[i,] = c(NewCor_GM,NewCor_K)
	}
	NumCorLess_GM = which(Correlation_GM < CorMat[,1])
	NumCorLess_K = which(Correlation_K < CorMat[,2])
	
	PVal_GM = length(NumCorLess_GM)/NumPermute
	PVal_K = length(NumCorLess_K)/NumPermute
	return(c(PVal_GM,PVal_K))
}

# GM_Retain.Ind = 1:54
# K_Retain.Ind = 1:54
###################################################################
# ContactNum.sd.threshold = 5
# (GM_Retain.Ind = which(GM_ContactNum.sts[,2] < ContactNum.sd.threshold))
# (K_Retain.Ind = which(K_ContactNum.sts[,2] < ContactNum.sd.threshold))
# 
GM_Retain.Ind = 1:54
K_Retain.Ind = 1:54
IgnoreNodes = c(1,11,16,19,20,28,30,31,39,44,45,51)
(GM_Retain.Ind = setdiff(GM_Retain.Ind,IgnoreNodes))
(K_Retain.Ind = setdiff(K_Retain.Ind, IgnoreNodes))

# GetPval(GM_Retain.Ind, DNase_Duke_GM, CI_GM, K_Retain.Ind, DNase_Duke_K, CI_K)
# GetPval(GM_Retain.Ind, CTCF_GM, CI_GM, K_Retain.Ind, CTCF_K, CI_K)
# GetPval(GM_Retain.Ind, H3K4me1_GM, CI_GM, K_Retain.Ind, H3K4me1_K, CI_K)
# GetPval(GM_Retain.Ind, H3K4me2_GM, CI_GM, K_Retain.Ind, H3K4me2_K, CI_K)
# GetPval(GM_Retain.Ind, H3K4me3_GM, CI_GM, K_Retain.Ind, H3K4me3_K, CI_K)
# GetPval(GM_Retain.Ind, H3K9ac_GM, CI_GM, K_Retain.Ind, H3K9ac_K, CI_K)
# GetPval(GM_Retain.Ind, H4K20me1_GM, CI_GM, K_Retain.Ind, H4K20me1_K, CI_K)
# GetPval(GM_Retain.Ind, H3K27me3_GM, CI_GM, K_Retain.Ind, H3K27me3_K, CI_K)
# GetPval(GM_Retain.Ind, H3K27ac_GM, CI_GM, K_Retain.Ind, H3K27ac_K, CI_K)
# GetPval(GM_Retain.Ind, H3K36me3_GM, CI_GM, K_Retain.Ind, H3K36me3_K, CI_K)

###################################################################
###################################################################
###################################################################


###################################################################


CI_GM_Normalize = CI_GM/max(CI_GM)
CI_K_Normalize = CI_K/max(CI_K)
DNase_Duke_GM_Normalize = DNase_Duke_GM/max(DNase_Duke_GM)
DNase_Duke_K_Normalize = DNase_Duke_K/max(DNase_Duke_K)
DNase_UW_GM_Normalize = DNase_UW_GM/max(DNase_UW_GM)
DNase_UW_K_Normalize = DNase_UW_K/max(DNase_UW_K)
H3K4me1_GM_Normalize = H3K4me1_GM/max(H3K4me1_GM)
H3K4me1_K_Normalize = H3K4me1_K/max(H3K4me1_K)
CTCF_GM_Normalize = CTCF_GM/max(CTCF_GM)
CTCF_K_Normalize = CTCF_K/max(CTCF_K)

# ----------------------------------------------------------------------------------------
H3K4me3_Broad_GM_Normalize = H3K4me3_Broad_GM/max(H3K4me3_Broad_GM)
H3K4me3_Broad_K_Normalize = H3K4me3_Broad_K/max(H3K4me3_Broad_K)
H3K4me3_UW_GM_Normalize = H3K4me3_UW_GM/max(H3K4me3_UW_GM)
H3K4me3_UW_K_Normalize = H3K4me3_UW_K/max(H3K4me3_UW_K)

CTCF_Broad_GM_Normalize = CTCF_Broad_GM/max(CTCF_Broad_GM)
CTCF_Broad_K_Normalize = CTCF_Broad_K/max(CTCF_Broad_K)
CTCF_UW_GM_Normalize = CTCF_UW_GM/max(CTCF_UW_GM)
CTCF_UW_K_Normalize = CTCF_UW_K/max(CTCF_UW_K)
CTCF_Uta_GM_Normalize = CTCF_Uta_GM/max(CTCF_Uta_GM)
CTCF_Uta_K_Normalize = CTCF_Uta_K/max(CTCF_Uta_K)
# ----------------------------------------------------------------------------------------



#----------
DNase_Duke_GM_nozeroind = which(DNase_Duke_GM_Normalize != 0)
DNase_Duke_GM_sign = GetSign(DNase_Duke_GM_Normalize[DNase_Duke_GM_nozeroind])
CI_GM_sign = GetSign(CI_GM_Normalize[DNase_Duke_GM_nozeroind])
rbind(DNase_Duke_GM_nozeroind, DNase_Duke_GM_sign - CI_GM_sign)
DNase_Duke_GM_sametrend = CheckSameTrend(DNase_Duke_GM_nozeroind,DNase_Duke_GM_sign - CI_GM_sign)
DNase_Duke_GM_Cor = GetRegionCorrelation(DNase_Duke_GM_Normalize,CI_GM_Normalize,DNase_Duke_GM_sametrend)

DNase_Duke_K_nozeroind = which(DNase_Duke_K_Normalize != 0)
DNase_Duke_K_sign = GetSign(DNase_Duke_K_Normalize[DNase_Duke_K_nozeroind])
CI_K_sign = GetSign(CI_K_Normalize[DNase_Duke_K_nozeroind])
rbind(DNase_Duke_K_nozeroind, DNase_Duke_K_sign - CI_K_sign)
DNase_Duke_K_sametrend = CheckSameTrend(DNase_Duke_K_nozeroind,DNase_Duke_K_sign - CI_K_sign)
DNase_Duke_K_Cor = GetRegionCorrelation(DNase_Duke_K_Normalize,CI_K_Normalize,DNase_Duke_K_sametrend)

#----------
DNase_UW_GM_nozeroind = which(DNase_UW_GM_Normalize != 0)
DNase_UW_GM_sign = GetSign(DNase_UW_GM_Normalize[DNase_UW_GM_nozeroind])
CI_GM_sign = GetSign(CI_GM_Normalize[DNase_UW_GM_nozeroind])
rbind(DNase_UW_GM_nozeroind, DNase_UW_GM_sign - CI_GM_sign)
DNase_UW_GM_sametrend = CheckSameTrend(DNase_UW_GM_nozeroind,DNase_UW_GM_sign - CI_GM_sign)
DNase_UW_GM_Cor = GetRegionCorrelation(DNase_UW_GM_Normalize,CI_GM_Normalize,DNase_UW_GM_sametrend)


DNase_UW_K_nozeroind = which(DNase_UW_K_Normalize != 0)
DNase_UW_K_sign = GetSign(DNase_UW_K_Normalize[DNase_UW_K_nozeroind])
CI_K_sign = GetSign(CI_K_Normalize[DNase_UW_K_nozeroind])
rbind(DNase_UW_K_nozeroind, DNase_UW_K_sign - CI_K_sign)
DNase_UW_K_sametrend = CheckSameTrend(DNase_UW_K_nozeroind,DNase_UW_K_sign - CI_K_sign)
DNase_UW_K_Cor = GetRegionCorrelation(DNase_UW_K_Normalize,CI_K_Normalize,DNase_UW_K_sametrend)

#----------
H3K4me1_GM_nozeroind = which(H3K4me1_GM_Normalize != 0)
H3K4me1_GM_sign = GetSign(H3K4me1_GM_Normalize[H3K4me1_GM_nozeroind])
CI_GM_sign = GetSign(CI_GM_Normalize[H3K4me1_GM_nozeroind])
rbind(H3K4me1_GM_nozeroind,H3K4me1_GM_sign - CI_GM_sign)
H3K4me1_GM_sametrend = CheckSameTrend(H3K4me1_GM_nozeroind,H3K4me1_GM_sign - CI_GM_sign)
H3K4me1_GM_Cor = GetRegionCorrelation(H3K4me1_GM_Normalize,CI_GM_Normalize,H3K4me1_GM_sametrend)

H3K4me1_K_nozeroind = which(H3K4me1_K_Normalize != 0)
H3K4me1_K_sign = GetSign(H3K4me1_K_Normalize[H3K4me1_K_nozeroind])
CI_K_sign = GetSign(CI_K_Normalize[H3K4me1_K_nozeroind])
rbind(H3K4me1_K_nozeroind,H3K4me1_K_sign - CI_K_sign)
H3K4me1_K_sametrend = CheckSameTrend(H3K4me1_K_nozeroind,H3K4me1_K_sign - CI_K_sign)
H3K4me1_K_Cor = GetRegionCorrelation(H3K4me1_K_Normalize,CI_K_Normalize,H3K4me1_K_sametrend)

#----------
CTCF_GM_nozeroind = which(CTCF_GM_Normalize != 0)
CTCF_GM_sign = GetSign(CTCF_GM_Normalize[CTCF_GM_nozeroind])
CI_GM_sign = GetSign(CI_GM_Normalize[CTCF_GM_nozeroind])
rbind(CTCF_GM_nozeroind, CTCF_GM_sign - CI_GM_sign)
CTCF_GM_sametrend = CheckSameTrend(CTCF_GM_nozeroind,CTCF_GM_sign - CI_GM_sign)
CTCF_GM_Cor = GetRegionCorrelation(CTCF_GM_Normalize,CI_GM_Normalize,CTCF_GM_sametrend)

CTCF_K_nozeroind = which(CTCF_K_Normalize != 0)
CTCF_K_sign = GetSign(CTCF_K_Normalize[CTCF_K_nozeroind])
CI_K_sign = GetSign(CI_K_Normalize[CTCF_K_nozeroind])
rbind(CTCF_K_nozeroind, CTCF_K_sign - CI_K_sign)
CTCF_K_sametrend = CheckSameTrend(CTCF_K_nozeroind,CTCF_K_sign - CI_K_sign)
CTCF_K_Cor = GetRegionCorrelation(CTCF_K_Normalize,CI_K_Normalize,CTCF_K_sametrend)

###################################################################
library(Matching)
ks.boot(DNase_Duke_GM_Normalize[DNase_Duke_GM_nozeroind], CI_GM_Normalize[DNase_Duke_GM_nozeroind])
ks.boot(DNase_Duke_K_Normalize[DNase_Duke_K_nozeroind], CI_K_Normalize[DNase_Duke_K_nozeroind])

ks.boot(DNase_UW_GM_Normalize[DNase_UW_GM_nozeroind], CI_GM_Normalize[DNase_UW_GM_nozeroind])
ks.boot(DNase_UW_K_Normalize[DNase_UW_K_nozeroind], CI_K_Normalize[DNase_UW_K_nozeroind])

ks.boot(H3K4me1_GM_Normalize[H3K4me1_GM_nozeroind], CI_GM_Normalize[H3K4me1_GM_nozeroind])
ks.boot(H3K4me1_K_Normalize[H3K4me1_K_nozeroind], CI_K_Normalize[H3K4me1_K_nozeroind])

###################################################################
###################################################################
###################################################################
###################################################################
###################################################################
# #------------------------------------------------------------------
library(reshape)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(gridExtra)
# cairo_ps(file="figs/Bar.ProbNodesSummation.ps",width=10,height=3)
# pdf(file="/tmp/Bar.Median.Neighbours.Number.pdf",width=20,height=6)
# # #------------------------------------------------------------------
# YMax = max(GM_ContactNum.sts[,1], K_ContactNum.sts[,1])
# df = data.frame(Node = 1:NumNodes, GM12878 = GM_ContactNum.sts[,1], K562 = K_ContactNum.sts[,1])
# df_long = melt(df,id="Node")
# MyPal = brewer.pal(9,"Set1")
# dodge <- position_dodge(width=0.5)
# gp1 <- ggplot(df_long, aes(x= factor(Node,levels=1:NumNodes), y=value, color= variable,group=variable))
# gp1 = gp1 +
#   theme_bw() +
#   theme(
#     text = element_text(size=15),
#     axis.text.x=element_text(angle=30, hjust=1,vjust=1),
#     legend.justification=c(1,1),
#     legend.position=c(1,1) ) +
#   labs(x= paste("Primer sites"),y = "Median of neighbour number") +
#   scale_x_discrete(label=xlabel) +
#   scale_color_manual("",
#     values = c( "GM12878" = MyPal[1], "K562" = MyPal[2]))
# gp1 = gp1 + geom_linerange(
# 	subset=.(variable %in% c("GM12878","K562")), 
# 	aes(x= factor(Node,levels=1:NumNodes),ymin=0,ymax=value,color=variable), position=  dodge,size=2)
# # gp + geom_bar(subset=.(variable %in% c("Random")), stat="identity",fill=I("grey70"),alpha=0.8)
# # gp
# # dev.off()
# 
# 
# # pdf(file="/tmp/Bar.Sd.Neighbours.Number.pdf",width=20,height=6)
# df = data.frame(Node = 1:NumNodes, GM12878 = GM_ContactNum.sts[,2], K562 = K_ContactNum.sts[,2])
# df_long = melt(df,id="Node")
# MyPal = brewer.pal(9,"Set1")
# dodge <- position_dodge(width=0.5)
# gp2 <- ggplot(df_long, aes(x= factor(Node,levels=1:NumNodes), y=value, color= variable,group=variable))
# gp2 = gp2 +
#   theme_bw() +
#   theme(
#     text = element_text(size=15),
#     axis.text.x=element_text(angle=30, hjust=1,vjust=1),
#     legend.justification=c(1,1),
#     legend.position=c(1,1) ) +
#   labs(x= paste("Primer sites"),y = "Sd of neighbour number") +
#   scale_x_discrete(label=xlabel) +
#   scale_color_manual("",
#     values = c( "GM12878" = MyPal[1], "K562" = MyPal[2]))
# gp2 = gp2 + geom_linerange(
# 	subset=.(variable %in% c("GM12878","K562")), 
# 	aes(x= factor(Node,levels=1:NumNodes),ymin=0,ymax=value,color=variable), position=  dodge,size=2)
# gp2 = gp2 + coord_cartesian(ylim=c(0,YMax))
# # gp + geom_bar(subset=.(variable %in% c("Random")), stat="identity",fill=I("grey70"),alpha=0.8)
# # gp
# # dev.off()
# 
# pdf(file="/tmp/Bar.Neighbours.Number.pdf",width=20,height=12)
# grid.arrange(gp1,gp2)
# dev.off()

#------------------------------------------------------------------------------
# DNase Duke
df_ignore =data.frame(Node=setdiff(1:NumNodes,DNase_Duke_GM_nozeroind), variable=factor("IgnoreNodes"),value = 0)

df = data.frame(Node = 1:NumNodes, DNase_Duke_GM = DNase_Duke_GM_Normalize, CI_GM = CI_GM_Normalize)
df_long = melt(df,id="Node")
df_long = rbind(df_long,df_ignore)

MyPal = brewer.pal(9,"Set1")
dodge <- position_dodge(width=0.5)
gp1 <- ggplot(df_long, aes(x= factor(Node,levels=1:NumNodes), y=value, color= variable,group=variable))
gp1 = gp1 +
  theme_bw() +
  theme(
    text = element_text(size=15),
    axis.text.x=element_text(angle=30, hjust=1,vjust=1),
    legend.justification=c(1,1),
    legend.position=c(1,.9) ) +
  labs(x= paste("Primer sites"),y = "Normalized enrichment") +
  scale_x_discrete(label=xlabel) +
  scale_color_manual("",
    values = c( "DNase_Duke_GM" = MyPal[1], "CI_GM" = MyPal[3], "IgnoreNodes" = "black"))
gp1 = gp1 + geom_linerange(
	subset=.(variable %in% c("DNase_Duke_GM","CI_GM")), 
	aes(x= factor(Node,levels=1:NumNodes),ymin=0,ymax=value,color=variable), position=  dodge,size=2)
gp1 = gp1 + geom_point(	data=subset(df_long,variable %in% c("IgnoreNodes")),shape = 2,size=5)

for (i in 1:nrow(DNase_Duke_GM_sametrend)){
	gp1 = gp1 + geom_rect(
		data=data.frame (xmin=DNase_Duke_GM_sametrend[i,1]-0.5, xmax=DNase_Duke_GM_sametrend[i,2]+0.5, ymin=-Inf, ymax=Inf),
		aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
		color="grey50", 
		alpha=0.1, 
		inherit.aes = FALSE)

	gp1 = gp1 + geom_text(
		data =
    data.frame(x=mean(DNase_Duke_GM_sametrend[i,]),y=1,label=paste("cor=",sprintf("%.2f",DNase_Duke_GM_Cor[i]))),
		aes(x=x,y=y,label=label),
		inherit.aes = FALSE	)

		RegionEnds = c(DNase_Duke_GM_sametrend[i,1]:DNase_Duke_GM_sametrend[i,2])
		RegionData = DNase_Duke_GM_Normalize[RegionEnds]
		RegionInd = which(RegionData != 0)
		RegionEnds = RegionEnds[RegionInd]
		RegionData = RegionData[RegionInd]
		RegionData_CI = CI_GM_Normalize[RegionEnds]
		gp1 = gp1 + geom_line(
			data=data.frame (x=RegionEnds,y=RegionData),
			aes(x=x,y=y), 
			color=MyPal[1], 
			inherit.aes = FALSE)
		gp1 = gp1 + geom_line(
			data=data.frame (x=RegionEnds,y=RegionData_CI),
			aes(x=x,y=y), 
			color=MyPal[3], 
			inherit.aes = FALSE)
}




df_ignore =data.frame(Node=setdiff(1:NumNodes,DNase_Duke_K_nozeroind), variable=factor("IgnoreNodes"),value = 0)
df = data.frame(Node = 1:NumNodes, DNase_Duke_K = DNase_Duke_K_Normalize, CI_K = CI_K_Normalize)
df_long = melt(df,id="Node")
df_long = rbind(df_long,df_ignore)

MyPal = brewer.pal(9,"Set1")
dodge <- position_dodge(width=0.5)
gp2 <- ggplot(df_long, aes(x= factor(Node,levels=1:NumNodes), y=value, color= variable,group=variable))
gp2 = gp2 +
  theme_bw() +
  theme(
    text = element_text(size=15),
    axis.text.x=element_text(angle=30, hjust=1,vjust=1),
    legend.justification=c(1,1),
    legend.position=c(1,.9) ) +
  labs(x= paste("Primer sites"),y = "Normalized enrichment") +
  scale_x_discrete(label=xlabel) +
  scale_color_manual("",
    values = c( "DNase_Duke_K" = MyPal[2], "CI_K" = MyPal[3], "IgnoreNodes" = "black"))
gp2 = gp2 + geom_linerange(
	subset=.(variable %in% c("DNase_Duke_K","CI_K")), 
	aes(x= factor(Node,levels=1:NumNodes),ymin=0,ymax=value,color=variable), position=  dodge,size=2)
gp2 = gp2 + geom_point(	data=subset(df_long,variable %in% c("IgnoreNodes")),shape = 2,size=5)

for (i in 1:nrow(DNase_Duke_K_sametrend)){
	gp2 = gp2 + geom_rect(
		data=data.frame (xmin=DNase_Duke_K_sametrend[i,1]-0.5, xmax=DNase_Duke_K_sametrend[i,2]+0.5, ymin=-Inf, ymax=Inf),
		aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
		color="grey50", 
		alpha=0.1, 
		inherit.aes = FALSE)

	gp2 = gp2 + geom_text(
		data =
    data.frame(x=mean(DNase_Duke_K_sametrend[i,]),y=1,label=paste("cor=",sprintf("%.2f",DNase_Duke_K_Cor[i]))),
		aes(x=x,y=y,label=label),
		inherit.aes = FALSE	)

	RegionEnds = c(DNase_Duke_K_sametrend[i,1]:DNase_Duke_K_sametrend[i,2])
	RegionData = DNase_Duke_K_Normalize[RegionEnds]
	RegionInd = which(RegionData != 0)
	RegionEnds = RegionEnds[RegionInd]
	RegionData = RegionData[RegionInd]
	RegionData_CI = CI_K_Normalize[RegionEnds]
	gp2 = gp2 + geom_line(
		data=data.frame (x=RegionEnds,y=RegionData),
		aes(x=x,y=y), 
		color=MyPal[2], 
		inherit.aes = FALSE)
	gp2 = gp2 + geom_line(
		data=data.frame (x=RegionEnds,y=RegionData_CI),
		aes(x=x,y=y), 
		color=MyPal[3], 
		inherit.aes = FALSE)
}


pdf(file="/tmp/Bar.DNase.CI.Duke.pdf",width=20,height=12)
grid.arrange(gp1,gp2)
dev.off()

#------------------------------------------------------------------------------
df_ignore =data.frame(Node=setdiff(1:NumNodes,DNase_UW_GM_nozeroind), variable=factor("IgnoreNodes"),value = 0)

df = data.frame(Node = 1:NumNodes, DNase_UW_GM = DNase_UW_GM_Normalize, CI_GM = CI_GM_Normalize)
df_long = melt(df,id="Node")
df_long = rbind(df_long,df_ignore)

MyPal = brewer.pal(9,"Set1")
dodge <- position_dodge(width=0.5)
gp1 <- ggplot(df_long, aes(x= factor(Node,levels=1:NumNodes), y=value, color= variable,group=variable))
gp1 = gp1 +
  theme_bw() +
  theme(
    text = element_text(size=15),
    axis.text.x=element_text(angle=30, hjust=1,vjust=1),
    legend.justification=c(1,1),
    legend.position=c(1,.9) ) +
  labs(x= paste("Primer sites"),y = "Normalized enrichment") +
  scale_x_discrete(label=xlabel) +
  scale_color_manual("",
    values = c( "DNase_UW_GM" = MyPal[1], "CI_GM" = MyPal[3], "IgnoreNodes" = "black"))
gp1 = gp1 + geom_linerange(
	subset=.(variable %in% c("DNase_UW_GM","CI_GM")), 
	aes(x= factor(Node,levels=1:NumNodes),ymin=0,ymax=value,color=variable), position=  dodge,size=2)
gp1 = gp1 + geom_point(	data=subset(df_long,variable %in% c("IgnoreNodes")),shape = 2,size=5)

for (i in 1:nrow(DNase_UW_GM_sametrend)){
	gp1 = gp1 + geom_rect(
		data=data.frame (xmin=DNase_UW_GM_sametrend[i,1]-0.5, xmax=DNase_UW_GM_sametrend[i,2]+0.5, ymin=-Inf, ymax=Inf),
		aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
		color="grey50", 
		alpha=0.1, 
		inherit.aes = FALSE)
	
	gp1 = gp1 + geom_text(
		data =
    data.frame(x=mean(DNase_UW_GM_sametrend[i,]),y=1,label=paste("cor=",sprintf("%.2f",DNase_UW_GM_Cor[i]))),
		aes(x=x,y=y,label=label),
		inherit.aes = FALSE	)

	RegionEnds = c(DNase_UW_GM_sametrend[i,1]:DNase_UW_GM_sametrend[i,2])
	RegionData = DNase_UW_GM_Normalize[RegionEnds]
	RegionInd = which(RegionData != 0)
	RegionEnds = RegionEnds[RegionInd]
	RegionData = RegionData[RegionInd]
	RegionData_CI = CI_GM_Normalize[RegionEnds]

	gp1 = gp1 + geom_line(
		data=data.frame (x=RegionEnds,y=RegionData),
		aes(x=x,y=y), 
		color=MyPal[1], 
		inherit.aes = FALSE)
	gp1 = gp1 + geom_line(
		data=data.frame (x=RegionEnds,y=RegionData_CI),
		aes(x=x,y=y), 
		color=MyPal[3], 
		inherit.aes = FALSE)

}



df_ignore =data.frame(Node=setdiff(1:NumNodes,DNase_UW_K_nozeroind), variable=factor("IgnoreNodes"),value = 0)
df = data.frame(Node = 1:NumNodes, DNase_UW_K = DNase_UW_K_Normalize, CI_K = CI_K_Normalize)
df_long = melt(df,id="Node")
df_long = rbind(df_long,df_ignore)

MyPal = brewer.pal(9,"Set1")
dodge <- position_dodge(width=0.5)
gp2 <- ggplot(df_long, aes(x= factor(Node,levels=1:NumNodes), y=value, color= variable,group=variable))
gp2 = gp2 +
  theme_bw() +
  theme(
    text = element_text(size=15),
    axis.text.x=element_text(angle=30, hjust=1,vjust=1),
    legend.justification=c(1,1),
    legend.position=c(1,.9) ) +
  labs(x= paste("Primer sites"),y = "Normalized enrichment") +
  scale_x_discrete(label=xlabel) +
  scale_color_manual("",
    values = c( "DNase_UW_K" = MyPal[2], "CI_K" = MyPal[3], "IgnoreNodes" = "black"))
gp2 = gp2 + geom_linerange(
	subset=.(variable %in% c("DNase_UW_K","CI_K")), 
	aes(x= factor(Node,levels=1:NumNodes),ymin=0,ymax=value,color=variable), position=  dodge,size=2)
gp2 = gp2 + geom_point(	data=subset(df_long,variable %in% c("IgnoreNodes")),shape = 2,size=5)

for (i in 1:nrow(DNase_UW_K_sametrend)){
	gp2 = gp2 + geom_rect(
		data=data.frame (xmin=DNase_UW_K_sametrend[i,1]-0.5, xmax=DNase_UW_K_sametrend[i,2]+0.5, ymin=-Inf, ymax=Inf),
		aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
		color="grey50", 
		alpha=0.1, 
		inherit.aes = FALSE)

	gp2 = gp2 + geom_text(
		data =
    data.frame(x=mean(DNase_UW_K_sametrend[i,]),y=1,label=paste("cor=",sprintf("%.2f",DNase_UW_K_Cor[i]))),
		aes(x=x,y=y,label=label),
		inherit.aes = FALSE	)

		RegionEnds = c(DNase_UW_K_sametrend[i,1]:DNase_UW_K_sametrend[i,2])
		RegionData = DNase_UW_K_Normalize[RegionEnds]
		RegionInd = which(RegionData != 0)
		RegionEnds = RegionEnds[RegionInd]
		RegionData = RegionData[RegionInd]
		RegionData_CI = CI_K_Normalize[RegionEnds]
		gp2 = gp2 + geom_line(
			data=data.frame (x=RegionEnds,y=RegionData),
			aes(x=x,y=y), 
			color=MyPal[2], 
			inherit.aes = FALSE)
		gp2 = gp2 + geom_line(
			data=data.frame (x=RegionEnds,y=RegionData_CI),
			aes(x=x,y=y), 
			color=MyPal[3], 
			inherit.aes = FALSE)		
}

pdf(file="/tmp/Bar.DNase.CI.UW.pdf",width=20,height=12)
grid.arrange(gp1,gp2)
dev.off()
#------------------------------------------------------------------------------



df_ignore =data.frame(Node=setdiff(1:NumNodes,H3K4me1_GM_nozeroind), variable=factor("IgnoreNodes"),value = 0)

df = data.frame(Node = 1:NumNodes, H3K4me1_GM = H3K4me1_GM_Normalize, CI_GM = CI_GM_Normalize)
df_long = melt(df,id="Node")
df_long = rbind(df_long,df_ignore)

MyPal = brewer.pal(9,"Set1")
dodge <- position_dodge(width=0.5)
gp1 <- ggplot(df_long, aes(x= factor(Node,levels=1:NumNodes), y=value, color= variable,group=variable))
gp1 = gp1 +
  theme_bw() +
  theme(
    text = element_text(size=15),
    axis.text.x=element_text(angle=30, hjust=1,vjust=1),
    legend.justification=c(1,1),
    legend.position=c(1,.9) ) +
  labs(x= paste("Primer sites"),y = "Normalized enrichment") +
  scale_x_discrete(label=xlabel) +
  scale_color_manual("",
    values = c( "H3K4me1_GM" = MyPal[1], "CI_GM" = MyPal[3], "IgnoreNodes" = "black"))
gp1 = gp1 + geom_linerange(
	subset=.(variable %in% c("H3K4me1_GM","CI_GM")), 
	aes(x= factor(Node,levels=1:NumNodes),ymin=0,ymax=value,color=variable), position=  dodge,size=2)
gp1 = gp1 + geom_point(	data=subset(df_long,variable %in% c("IgnoreNodes")),shape = 2,size=5)

for (i in 1:nrow(H3K4me1_GM_sametrend)){
	gp1 = gp1 + geom_rect(
		data=data.frame (xmin=H3K4me1_GM_sametrend[i,1]-0.5, xmax=H3K4me1_GM_sametrend[i,2]+0.5, ymin=-Inf, ymax=Inf),
		aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
		color="grey50", 
		alpha=0.1, 
		inherit.aes = FALSE)

	gp1 = gp1 + geom_text(
		data =
    data.frame(x=mean(H3K4me1_GM_sametrend[i,]),y=1,label=paste("cor=",sprintf("%.2f",H3K4me1_GM_Cor[i]))),
		aes(x=x,y=y,label=label),
		inherit.aes = FALSE	)

		RegionEnds = c(H3K4me1_GM_sametrend[i,1]:H3K4me1_GM_sametrend[i,2])
		RegionData = H3K4me1_GM_Normalize[RegionEnds]
		RegionInd = which(RegionData != 0)
		RegionEnds = RegionEnds[RegionInd]
		RegionData = RegionData[RegionInd]
		RegionData_CI = CI_GM_Normalize[RegionEnds]
		gp1 = gp1 + geom_line(
			data=data.frame (x=RegionEnds,y=RegionData),
			aes(x=x,y=y), 
			color=MyPal[1], 
			inherit.aes = FALSE)
		gp1 = gp1 + geom_line(
			data=data.frame (x=RegionEnds,y=RegionData_CI),
			aes(x=x,y=y), 
			color=MyPal[3], 
			inherit.aes = FALSE)
}



df_ignore =data.frame(Node=setdiff(1:NumNodes,H3K4me1_K_nozeroind), variable=factor("IgnoreNodes"),value = 0)
df = data.frame(Node = 1:NumNodes, H3K4me1_K = H3K4me1_K/max(H3K4me1_K), CI_K = CI_K/max(CI_K))
df_long = melt(df,id="Node")
df_long = rbind(df_long,df_ignore)

MyPal = brewer.pal(9,"Set1")
dodge <- position_dodge(width=0.5)
gp2 <- ggplot(df_long, aes(x= factor(Node,levels=1:NumNodes), y=value, color= variable,group=variable))
gp2 = gp2 +
  theme_bw() +
  theme(
    text = element_text(size=15),
    axis.text.x=element_text(angle=30, hjust=1,vjust=1),
    legend.justification=c(1,1),
    legend.position=c(1,.9) ) +
  labs(x= paste("Primer sites"),y = "Normalized enrichment") +
  scale_x_discrete(label=xlabel) +
  scale_color_manual("",
    values = c( "H3K4me1_K" = MyPal[2], "CI_K" = MyPal[3], "IgnoreNodes" = "black"))
gp2 = gp2 + geom_linerange(
	subset=.(variable %in% c("H3K4me1_K","CI_K")), 
	aes(x= factor(Node,levels=1:NumNodes),ymin=0,ymax=value,color=variable), position=  dodge,size=2)
gp2 = gp2 + geom_point(	data=subset(df_long,variable %in% c("IgnoreNodes")),shape = 2,size=5)
for (i in 1:nrow(H3K4me1_K_sametrend)){
	gp2 = gp2 + geom_rect(
		data=data.frame (xmin=H3K4me1_K_sametrend[i,1]-0.5, xmax=H3K4me1_K_sametrend[i,2]+0.5, ymin=-Inf, ymax=Inf),
		aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
		color="grey50", 
		alpha=0.1, 
		inherit.aes = FALSE)

	gp2 = gp2 + geom_text(
		data =
    data.frame(x=mean(H3K4me1_K_sametrend[i,]),y=1,label=paste("cor=",sprintf("%.2f",H3K4me1_K_Cor[i]))),
		aes(x=x,y=y,label=label),
		inherit.aes = FALSE	)

		RegionEnds = c(H3K4me1_K_sametrend[i,1]:H3K4me1_K_sametrend[i,2])
		RegionData = H3K4me1_K_Normalize[RegionEnds]
		RegionInd = which(RegionData != 0)
		RegionEnds = RegionEnds[RegionInd]
		RegionData = RegionData[RegionInd]
		RegionData_CI = CI_K_Normalize[RegionEnds]
		gp2 = gp2 + geom_line(
			data=data.frame (x=RegionEnds,y=RegionData),
			aes(x=x,y=y), 
			color=MyPal[2], 
			inherit.aes = FALSE)
		gp2 = gp2 + geom_line(
			data=data.frame (x=RegionEnds,y=RegionData_CI),
			aes(x=x,y=y), 
			color=MyPal[3], 
			inherit.aes = FALSE)
}


pdf(file="/tmp/Bar.H3K4me1.CI.pdf",width=20,height=12)
grid.arrange(gp1,gp2)
dev.off()




#------------------------------------------------------------------------------
df_ignore =data.frame(Node=setdiff(1:NumNodes,CTCF_GM_nozeroind), variable=factor("IgnoreNodes"),value = 0)

df = data.frame(Node = 1:NumNodes, CTCF_GM = CTCF_GM_Normalize, CI_GM = CI_GM_Normalize)
df_long = melt(df,id="Node")
df_long = rbind(df_long,df_ignore)

MyPal = brewer.pal(9,"Set1")
dodge <- position_dodge(width=0.5)
gp1 <- ggplot(df_long, aes(x= factor(Node,levels=1:NumNodes), y=value, color= variable,group=variable))
gp1 = gp1 +
  theme_bw() +
  theme(
    text = element_text(size=15),
    axis.text.x=element_text(angle=30, hjust=1,vjust=1),
    legend.justification=c(1,1),
    legend.position=c(1,.9) ) +
  labs(x= paste("Primer sites"),y = "Normalized enrichment") +
  scale_x_discrete(label=xlabel) +
  scale_color_manual("",
    values = c( "CTCF_GM" = MyPal[1], "CI_GM" = MyPal[3], "IgnoreNodes" = "black"))
gp1 = gp1 + geom_linerange(
	subset=.(variable %in% c("CTCF_GM","CI_GM")), 
	aes(x= factor(Node,levels=1:NumNodes),ymin=0,ymax=value,color=variable), position=  dodge,size=2)
gp1 = gp1 + geom_point(	data=subset(df_long,variable %in% c("IgnoreNodes")),shape = 2,size=5)

for (i in 1:nrow(CTCF_GM_sametrend)){
	gp1 = gp1 + geom_rect(
		data=data.frame (xmin=CTCF_GM_sametrend[i,1]-0.5, xmax=CTCF_GM_sametrend[i,2]+0.5, ymin=-Inf, ymax=Inf),
		aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
		color="grey50", 
		alpha=0.1, 
		inherit.aes = FALSE)
	
	gp1 = gp1 + geom_text(
		data =
    data.frame(x=mean(CTCF_GM_sametrend[i,]),y=1,label=paste("cor=",sprintf("%.2f",CTCF_GM_Cor[i]))),
		aes(x=x,y=y,label=label),
		inherit.aes = FALSE	)

	RegionEnds = c(CTCF_GM_sametrend[i,1]:CTCF_GM_sametrend[i,2])
	RegionData = CTCF_GM_Normalize[RegionEnds]
	RegionInd = which(RegionData != 0)
	RegionEnds = RegionEnds[RegionInd]
	RegionData = RegionData[RegionInd]
	RegionData_CI = CI_GM_Normalize[RegionEnds]

	gp1 = gp1 + geom_line(
		data=data.frame (x=RegionEnds,y=RegionData),
		aes(x=x,y=y), 
		color=MyPal[1], 
		inherit.aes = FALSE)
	gp1 = gp1 + geom_line(
		data=data.frame (x=RegionEnds,y=RegionData_CI),
		aes(x=x,y=y), 
		color=MyPal[3], 
		inherit.aes = FALSE)

}



df_ignore =data.frame(Node=setdiff(1:NumNodes,CTCF_K_nozeroind), variable=factor("IgnoreNodes"),value = 0)
df = data.frame(Node = 1:NumNodes, CTCF_K = CTCF_K_Normalize, CI_K = CI_K_Normalize)
df_long = melt(df,id="Node")
df_long = rbind(df_long,df_ignore)

MyPal = brewer.pal(9,"Set1")
dodge <- position_dodge(width=0.5)
gp2 <- ggplot(df_long, aes(x= factor(Node,levels=1:NumNodes), y=value, color= variable,group=variable))
gp2 = gp2 +
  theme_bw() +
  theme(
    text = element_text(size=15),
    axis.text.x=element_text(angle=30, hjust=1,vjust=1),
    legend.justification=c(1,1),
    legend.position=c(1,.9) ) +
  labs(x= paste("Primer sites"),y = "Normalized enrichment") +
  scale_x_discrete(label=xlabel) +
  scale_color_manual("",
    values = c( "CTCF_K" = MyPal[2], "CI_K" = MyPal[3], "IgnoreNodes" = "black"))
gp2 = gp2 + geom_linerange(
	subset=.(variable %in% c("CTCF_K","CI_K")), 
	aes(x= factor(Node,levels=1:NumNodes),ymin=0,ymax=value,color=variable), position=  dodge,size=2)
gp2 = gp2 + geom_point(	data=subset(df_long,variable %in% c("IgnoreNodes")),shape = 2,size=5)

for (i in 1:nrow(CTCF_K_sametrend)){
	gp2 = gp2 + geom_rect(
		data=data.frame (xmin=CTCF_K_sametrend[i,1]-0.5, xmax=CTCF_K_sametrend[i,2]+0.5, ymin=-Inf, ymax=Inf),
		aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
		color="grey50", 
		alpha=0.1, 
		inherit.aes = FALSE)

	gp2 = gp2 + geom_text(
		data =
    data.frame(x=mean(CTCF_K_sametrend[i,]),y=1,label=paste("cor=",sprintf("%.2f",CTCF_K_Cor[i]))),
		aes(x=x,y=y,label=label),
		inherit.aes = FALSE	)

		RegionEnds = c(CTCF_K_sametrend[i,1]:CTCF_K_sametrend[i,2])
		RegionData = CTCF_K_Normalize[RegionEnds]
		RegionInd = which(RegionData != 0)
		RegionEnds = RegionEnds[RegionInd]
		RegionData = RegionData[RegionInd]
		RegionData_CI = CI_K_Normalize[RegionEnds]
		gp2 = gp2 + geom_line(
			data=data.frame (x=RegionEnds,y=RegionData),
			aes(x=x,y=y), 
			color=MyPal[2], 
			inherit.aes = FALSE)
		gp2 = gp2 + geom_line(
			data=data.frame (x=RegionEnds,y=RegionData_CI),
			aes(x=x,y=y), 
			color=MyPal[3], 
			inherit.aes = FALSE)		
}

pdf(file="/tmp/Bar.CTCF.CI.pdf",width=20,height=12)
grid.arrange(gp1,gp2)
dev.off()




#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------

# CTCFPij_GM = Pij_GM[which(Pij_GM[,1] %in% CTCF_GM_nozeroind & Pij_GM[,2] %in% CTCF_GM_nozeroind),]
# CTCFPij_K = Pij_K[which(Pij_K[,1] %in% CTCF_K_nozeroind & Pij_K[,2] %in% CTCF_K_nozeroind),]
# 
# CTCFPij_GM_Ind = unique(sort(c(CTCFPij_GM[,1],CTCFPij_GM[,2])))
# CTCFPij_GM_SumVec = matrix(0,nrow=length(CTCFPij_GM_Ind),ncol=2)
# for (i in 1:length(CTCFPij_GM_Ind)){
# 	CTCFPij_GM_SumVec[i,] = c(CTCFPij_GM_Ind[i],sum(CTCFPij_GM[which(CTCFPij_GM[,1] %in% CTCFPij_GM_Ind[i] | CTCFPij_GM[,2] %in% CTCFPij_GM_Ind[i]),3]))
# }
# 
# CTCFPij_GM_SumVec_Normalize = cbind(CTCFPij_GM_SumVec[,1],CTCFPij_GM_SumVec[,2]/max(CTCFPij_GM_SumVec[,2]))
# plot(CTCFPij_GM_SumVec_Normalize,type="o")
# lines(cbind(CTCFPij_GM_Ind,CTCF_GM_Normalize[CTCFPij_GM_Ind]),col="red")
# 
# 
# CTCFPij_K_Ind = unique(sort(c(CTCFPij_K[,1],CTCFPij_K[,2])))
# CTCFPij_K_SumVec = matrix(0,nrow=length(CTCFPij_K_Ind),ncol=2)
# for (i in 1:length(CTCFPij_K_Ind)){
# 	CTCFPij_K_SumVec[i,] = c(CTCFPij_K_Ind[i],sum(CTCFPij_K[which(CTCFPij_K[,1] %in% CTCFPij_K_Ind[i] | CTCFPij_K[,2] %in% CTCFPij_K_Ind[i]),3]))
# }
# 
# CTCFPij_K_SumVec_Normalize = cbind(CTCFPij_K_SumVec[,1],CTCFPij_K_SumVec[,2]/max(CTCFPij_K_SumVec[,2]))
# plot(CTCFPij_K_SumVec_Normalize,type="o")
# lines(cbind(CTCFPij_K_Ind,CTCF_K_Normalize[CTCFPij_K_Ind]),col="red")
# 

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# For CTCF
# CTCF_threshold = 20
# CTCF_GM_thresholdind = CTCF_GM_nozeroind[which(CTCF_GM[CTCF_GM_nozeroind]>= CTCF_threshold)]
# CTCF_K_thresholdind = CTCF_K_nozeroind[which(CTCF_K[CTCF_K_nozeroind]>= CTCF_threshold)]
# 
# # Pij_CTCF_GM_All = NULL
# # Exp_CTCF_GM_All = NULL
# for (i in 1:length(CTCF_GM_thresholdind)){
# 	iNode = CTCF_GM_thresholdind[i]
# 	iNode = 11
# 	CTCFNodesInd = which((Pij_GM[,1]==iNode & Pij_GM[,2] %in% CTCF_GM_thresholdind) | (Pij_GM[,2]==iNode & Pij_GM[,1] %in% CTCF_GM_thresholdind))
# 	Pij_CTCF_Node = Pij_GM[CTCFNodesInd,]
# 	Exp_CTCF_Node = matrix(0,nrow=nrow(Pij_CTCF_Node),ncol=ncol(Pij_CTCF_Node))
# 
# 	CTCF_Cor_NodeInd = NULL
# 	for ( j in 1:nrow(Pij_CTCF_Node)){
# 		Node1 = Pij_CTCF_Node[j,1]
# 		Node2 = Pij_CTCF_Node[j,2]
# 		if (Node1 != iNode){
# 			CTCF_Cor_NodeInd = c(CTCF_Cor_NodeInd, Node1)
# 		} else {
# 			CTCF_Cor_NodeInd = c(CTCF_Cor_NodeInd, Node2)
# 		}
# 		Exp_CTCF_Node[j,] = c(Node1,Node2,CTCF_GM[Node1]*CTCF_GM[Node2])
# 	}
# 	print (Pij_CTCF_Node)
# 	# Pij_CTCF_K_All = rbind(Pij_CTCF_K_All,Pij_CTCF_Node)
# 	# Exp_CTCF_K_All = rbind(Exp_CTCF_K_All,Exp_CTCF_Node)
# 	cat(iNode, cor(Exp_CTCF_Node[,3],Pij_CTCF_Node[,3]), "\n")
# 	
# 
# 	df = data.frame(
# 		Node = CTCF_Cor_NodeInd, 
# 		CTCF_GM = Exp_CTCF_Node[,3]/max(Exp_CTCF_Node[,3]), 
# 		Pij_GM = Pij_CTCF_Node[,3]/max(Pij_CTCF_Node[,3]))
# 	df_long = melt(df,id="Node")
# 
# 	MyPal = brewer.pal(9,"Set1")
# 	dodge <- position_dodge(width=0)
# 	gp1 <- ggplot(df_long, aes(x= Node, y=value, color= variable,group=variable))
# 	gp1 = gp1 +
# 	  theme_bw() +
# 	  theme(
# 	    text = element_text(size=15),
# 	    axis.text.x=element_text(angle=30, hjust=1,vjust=1),
# 	    legend.justification=c(1,1),
# 	    legend.position=c(1,.9) ) +
# 	  labs(x= paste("Primer sites (", xlabel[iNode], ")"),y = "Normalized enrichment") +
# 	  scale_x_continuous(breaks=1:NumNodes,labels=xlabel,limits=c(1,NumNodes)) +
# 		# coord_cartesian(xlim=c(0, NumNodes)) + 
# 	  scale_color_manual("",
# 	    values = c( "CTCF_GM" = MyPal[1], "Pij_GM" = MyPal[3]))
# 	gp1 = gp1 + geom_line(
# 		subset=.(variable %in% c("CTCF_GM","Pij_GM")), 
# 		aes(x= Node,ymin=0,ymax=value,color=variable), position=  dodge,size=2)
# 	gp1 = gp1 + geom_point(
# 		subset=.(variable %in% c("CTCF_GM","Pij_GM")), 
# 		aes(x= Node,ymin=0,ymax=value,color=variable), position=  dodge,size=4)
# 	
# }
# # cat("all",cor(Pij_CTCF_GM_All[,3],Exp_CTCF_GM_All[,3]),"\n")
# 
# 
# 
# # Pij_CTCF_K_All = NULL
# # Exp_CTCF_K_All = NULL
# for (i in 1:length(CTCF_K_thresholdind)){
# 	iNode = CTCF_K_thresholdind[i]
# 
# 	iNode = 40
# 	iNode = 52
# 	iNode = 43
# 	CTCFNodesInd = which((Pij_K[,1]==iNode & Pij_K[,2] %in% CTCF_K_thresholdind) | (Pij_K[,2]==iNode & Pij_K[,1] %in% CTCF_K_thresholdind))
# 	Pij_CTCF_Node = Pij_K[CTCFNodesInd,]
# 	Exp_CTCF_Node = matrix(0,nrow=nrow(Pij_CTCF_Node),ncol=ncol(Pij_CTCF_Node))
# 
# 	CTCF_Cor_NodeInd = NULL
# 	for ( j in 1:nrow(Pij_CTCF_Node)){
# 		Node1 = Pij_CTCF_Node[j,1]
# 		Node2 = Pij_CTCF_Node[j,2]
# 		if (Node1 != iNode){
# 			CTCF_Cor_NodeInd = c(CTCF_Cor_NodeInd, Node1)
# 		} else {
# 			CTCF_Cor_NodeInd = c(CTCF_Cor_NodeInd, Node2)
# 		}
# 		Exp_CTCF_Node[j,] = c(Node1,Node2,CTCF_K[Node1]*CTCF_K[Node2])
# 	}
# 	print(Pij_CTCF_Node)
# 	# Pij_CTCF_K_All = rbind(Pij_CTCF_K_All,Pij_CTCF_Node)
# 	# Exp_CTCF_K_All = rbind(Exp_CTCF_K_All,Exp_CTCF_Node)
# 	cat(iNode, cor(Exp_CTCF_Node[,3],Pij_CTCF_Node[,3]), "\n")
# 	
# 
# 	df = data.frame(
# 		Node = CTCF_Cor_NodeInd, 
# 		CTCF_K = Exp_CTCF_Node[,3]/max(Exp_CTCF_Node[,3]), 
# 		Pij_K = Pij_CTCF_Node[,3]/max(Pij_CTCF_Node[,3]))
# 	df_long = melt(df,id="Node")
# 
# 
# 	MyPal = brewer.pal(9,"Set1")
# 	dodge <- position_dodge(width=0)
# 	gp1 <- ggplot(df_long, aes(x= Node, y=value, color= variable,group=variable))
# 	gp1 = gp1 +
# 	  theme_bw() +
# 	  theme(
# 	    text = element_text(size=15),
# 	    axis.text.x=element_text(angle=30, hjust=1,vjust=1),
# 	    legend.justification=c(1,1),
# 	    legend.position=c(1,.9) ) +
# 	  labs(x= paste("Primer sites (", xlabel[iNode], ")"),y = "Normalized enrichment") +
# 	  scale_x_continuous(breaks=1:NumNodes,labels=xlabel,limits=c(1,NumNodes)) +
# 		# coord_cartesian(xlim=c(0, NumNodes)) + 
# 	  scale_color_manual("",
# 	    values = c( "CTCF_K" = MyPal[2], "Pij_K" = MyPal[3]))
# 	gp1 = gp1 + geom_line(
# 		subset=.(variable %in% c("CTCF_K","Pij_K")), 
# 		aes(x= Node,ymin=0,ymax=value,color=variable), position=  dodge,size=2)
# 	gp1 = gp1 + geom_point(
# 		subset=.(variable %in% c("CTCF_K","Pij_K")), 
# 		aes(x= Node,ymin=0,ymax=value,color=variable), position=  dodge,size=4)
# 
# 	
# }
# # cat("all",cor(Pij_CTCF_K_All[,3],Exp_CTCF_K_All[,3]),"\n")
# 
# 
# plot(Pij_CTCF_Node[,3]/max(Pij_CTCF_Node[,3]),type="o")
# lines(Exp_CTCF_Node[,3]/max(Exp_CTCF_Node[,3]),col="red")
# 


