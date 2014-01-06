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

of_name = paste(DestDir,"ENm008.","p",PersistenceLength,".equ.seg.len.txt",sep="")
write.table(file=of_name,sprintf("%.3f",Seg_Len),row.names=FALSE,col.names=FALSE,quote=FALSE)


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


#############################################################################
# Output segment distance, specific point
of_name = paste(DestDir,"ENm008.","p",PersistenceLength,".equ.contact.point.ind.txt",sep="")
write.table(file=of_name,Equ_Contact_Point_Ind,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
of_name = paste(DestDir,"GM12878.","p",PersistenceLength,".pos.count.txt",sep="")
write.table(file=of_name,CFData_1,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
of_name = paste(DestDir,"K562.","p",PersistenceLength,".pos.count.txt",sep="")
write.table(file=of_name,CFData_2,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
of_name = paste(DestDir,"GM12878.","p",PersistenceLength,".point.ind.count.txt",sep="")
write.table(file=of_name,CF_1.Vec,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
of_name = paste(DestDir,"K562.","p",PersistenceLength,".point.ind.count.txt",sep="")
write.table(file=of_name,CF_2.Vec,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
of_name = paste(DestDir,"ENm008.","p",PersistenceLength,".point.ind.pos.txt",sep="")
write.table(file=of_name,cbind(Points_Map_Ind,Points_Pos),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
#----------------------------------------------------------------------------
# q()
#############################################################################

GetClosestEquPointIndFromExpPointPos <- function(ExpPointPos){
	ExpPointInd = NULL
	for (i in 1:length(ExpPointPos)){
		Points_Diff = abs(Points_Pos - ExpPointPos[i])
		ExpPointInd = c(ExpPointInd,which(Points_Diff == min(Points_Diff)))
	}
	return (GetEquPointIndFromExpPointInd(ExpPointInd))
}

#############################################################################
#----------------------------------------------------------------------------

(    HS48s_NodeInd   = GetEquPointIndFromExpPointPos(91497)     )   # HS48 start
(    HS48e_NodeInd   = GetEquPointIndFromExpPointPos(100530)    )   # HS48 end
(    HS40s_NodeInd   = GetEquPointIndFromExpPointPos(100530)    )   # HS40 start
(    HS40e_NodeInd   = GetEquPointIndFromExpPointPos(104687)    )   # HS40 end
(    HS33s_NodeInd   = GetEquPointIndFromExpPointPos(109838)    )   # HS33 start
(    HS33e_NodeInd   = GetEquPointIndFromExpPointPos(120933)    )   # HS33 end
(    HS10s_NodeInd   = GetEquPointIndFromExpPointPos(131220)    )   # HS10 start
(    HS10e_NodeInd   = GetEquPointIndFromExpPointPos(134334)    )   # HS10 end
(    A2s_NodeInd     = GetEquPointIndFromExpPointPos(147782)    )   # A2 start
(    A2e_NodeInd     = GetEquPointIndFromExpPointPos(167103)    )   # A2 end

#All_Points_Pos[setdiff(c(1:length(All_Points_Pos)),GetEquPointIndFromExpPointPos(Points_Pos))]
cat(
"SampleSize        \t",SampleSize,          "\n"   ,
"HS48 start index  \t",HS48s_NodeInd        ,"\n"  ,
"HS48 end index    \t",HS48e_NodeInd        ,"\n"  ,
"HS40 start index  \t",HS40s_NodeInd        ,"\n"  ,
"HS40 end index    \t",HS40e_NodeInd        ,"\n"  ,
"HS33 start index  \t",HS33s_NodeInd        ,"\n"  ,
"HS33 end index    \t",HS33e_NodeInd        ,"\n"  ,
"HS10 start index  \t",HS10s_NodeInd        ,"\n"  ,
"HS10 end index    \t",HS10e_NodeInd        ,"\n"  ,
"A2 start index    \t",A2s_NodeInd          ,"\n"  ,
"A2 end index      \t",A2e_NodeInd          ,"\n"  ,
sep= ""
)
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
# Scale Contact frequency
#----------------------------------------------------------------------------
Con.Vec            = CF_1.Vec[,1:2]

CF_1.Vec.Scale     = CF_1.Vec
CF_1.Sumcount      = sum(CF_1.Vec[,3])
CF_1.Vec.Scale[,3] = CF_1.Vec.Scale[,3]/CF_1.Sumcount

CF_2.Vec.Scale     = CF_2.Vec
CF_2.Sumcount      = sum(CF_2.Vec[,3])
CF_2.Vec.Scale[,3] = CF_2.Vec.Scale[,3]/CF_2.Sumcount

#------------------------------------------------------------------
# Read reference data, and merge to experimental segment
#------------------------------------------------------------------

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

Contact.RF.Count.Vec =	apply(Contact.RF.Count.Vec,2,function(x)
			{ matchind = which(x==0)
				if(length(matchind)){x[matchind]= min(x[which(x!=0)])}
				return(x) })

RF.Vec = cbind(Con.Vec, Contact.RF.Count.Vec)
RF.Vec.Scale = cbind(Con.Vec, apply(Contact.RF.Count.Vec,2,function(x){x/sum(x)}))

# #------------------------------------------------------------------
# #------------------------------------------------------------------


#------------------------------------------------------------------
# Ensemble contact for circos
#------------------------------------------------------------------

#------------------------------------------------------------------
# Propensity
#------------------------------------------------------------------
PS_1.Vec = cbind(Con.Vec,
	apply(RF.Vec.Scale[,3:ncol(RF.Vec.Scale)], 2, function(x){CF_1.Vec.Scale[,3]/x}))
PS_2.Vec = cbind(Con.Vec,
	apply(RF.Vec.Scale[,3:ncol(RF.Vec.Scale)], 2, function(x){CF_2.Vec.Scale[,3]/x}))


#----------------------------------------
# circos of Propensity
#----------------------------------------

# iCount = 1
# LinkStr = ""
# for (i in 1:nrow(ConInd_Vec)){
#     iRow = ConInd_Vec[i,1]
#     iCol = ConInd_Vec[i,2]
#     if(PS.Vec[i,3]<0.001){next}
#     PositionStr1.1 = IntervalPoints[iRow]
#     PositionStr1.2 = IntervalPoints[iRow+1]
#     PositionStr2.1 = IntervalPoints[iCol]
#     PositionStr2.2 = IntervalPoints[iCol+1]
#     Field1 = paste("segdup",sprintf("%04d hs16", iCount),sep="")
#     freq = sprintf("thickness=%.3f",PS.Vec[i,3])
#     maxfreq = sprintf("id=%.2f",PS.Max,sep="")
#     params = paste(freq,",",maxfreq,sep="")
#     LinkStr = paste(LinkStr, paste(Field1," ",PositionStr2.1," ",PositionStr2.2," ",params,sep=""),sep="\n")
#     LinkStr = paste(LinkStr, paste(Field1," ",PositionStr1.1," ",PositionStr1.2," ",params,sep=""),sep="\n")
#     iCount = iCount + 1
# }
# 
# of_circos = paste(DestDir,CELL,".",SampleSize,".cut",Cutoff,".prop.txt",sep="")
# write.table(file=of_circos,LinkStr,quote=FALSE,row.names=FALSE,col.names=FALSE)

#------------------------------------------------------------------
#----------------------------------------
# histogram of HS40 in GM and K cell
#----------------------------------------
# PS.Spec.Ind = which(PS.Vec[,1]==A2SegInd | PS.Vec[,2]==A2SegInd)
# LinkStr = ""
# for (i in 1:(length(PS.Spec.Ind))){
#   if (PS.Vec[PS.Spec.Ind[i],1] == A2SegInd){
#     iRow = PS.Vec[PS.Spec.Ind[i],2]
#   }else{
#     iRow = PS.Vec[PS.Spec.Ind[i],1]
#   }
#   PositionStr1.1 = IntervalPoints[iRow]
#   PositionStr1.2 = IntervalPoints[iRow+1]
#   Field1 = "hs16"
#   params = sprintf("%.2f",PS.Vec[PS.Spec.Ind[i],3])
#   LinkStr = paste(LinkStr,paste(Field1," ",PositionStr1.1," ",PositionStr1.2," ",params,sep=""),sep="\n")
# }
# 
# of_circos = paste(DestDir,CELL,".",SampleSize,".cut",Cutoff,".a",alpha,".l",LengthDiscard,".PS.A2SegInd.txt",sep="")
# write.table(file=of_circos,LinkStr,quote=FALSE,row.names=FALSE,col.names=FALSE)
#----------------------------------------

# PS.Spec.Ind = which(PS.Vec[,1]==HS40SegInd | PS.Vec[,2]==HS40SegInd)
# LinkStr = ""
# for (i in 1:(length(PS.Spec.Ind))){
#   if (PS.Vec[PS.Spec.Ind[i],1] == HS40SegInd){
#     iRow = PS.Vec[PS.Spec.Ind[i],2]
#   }else{
#     iRow = PS.Vec[PS.Spec.Ind[i],1]
#   }
#   PositionStr1.1 = IntervalPoints[iRow]
#   PositionStr1.2 = IntervalPoints[iRow+1]
#   Field1 = "hs16"
#   params = sprintf("%.2f",PS.Vec[PS.Spec.Ind[i],3])
#   LinkStr = paste(LinkStr,paste(Field1," ",PositionStr1.1," ",PositionStr1.2," ",params,sep=""),sep="\n")
# }
# 
# of_circos = paste(DestDir,CELL,".",SampleSize,".cut",Cutoff,".a",alpha,".PS.HS40Seg.txt",sep="")
# write.table(file=of_circos,LinkStr,quote=FALSE,row.names=FALSE,col.names=FALSE)


#------------------------------------------------------------------
# Bootstrap for p-value
#------------------------------------------------------------------
AlphaInd.Vec = cbind(c(5,seq(10,100,by=10)), c(1:11)+2)
DistanceInd.Vec = cbind(c(840,1140),3:4)

CutOff = 840
CutColInd = DistanceInd.Vec[which(DistanceInd.Vec[,1]==CutOff),2]
CutOffPath = paste(dirname(getwd()),"/result/analysis/ENm008/cut/",SampleSize,"/c",CutOff,sep="")
ConFiles = dir(CutOffPath,pattern="*.con",recursive=TRUE,full.names=TRUE)

#ConFiles = ConFiles[1:10]
NumFiles = length(ConFiles)

Weight.Vec = rep(0,NumFiles)
ConMatrix = matrix(0,nrow=nrow(Con.Vec),ncol=NumFiles)
for (i in 1:NumFiles){
	ConFile = ConFiles[i]
	Weight.Vec[i] = read.table(ConFile,comment.char="",nrow=1)[,3]
	ConMatrix[,i] = read.table(ConFile)[,3]
}

BootstrapTimes = 1000
BootSampleSize = NumFiles
Bootstrap_1.Con = matrix(FALSE,nrow=nrow(Con.Vec),ncol=BootstrapTimes)
Bootstrap_2.Con = matrix(FALSE,nrow=nrow(Con.Vec),ncol=BootstrapTimes)

ptm <- proc.time()
for (i in 1:BootstrapTimes){
	Con.Samples.Ind = sample(1:NumFiles,BootSampleSize,replace=T)
	Con.Nominator = apply(ConMatrix[,Con.Samples.Ind],1,sum)
	Con.Denominator = sum(Con.Nominator)
	Con.q = Con.Nominator / Con.Denominator
	Bootstrap_1.Con[,i] = (CF_1.Vec.Scale[,3] < Con.q)
	Bootstrap_2.Con[,i] = (CF_2.Vec.Scale[,3] < Con.q)
	# print(i)
}
PVal_1.sim = apply(Bootstrap_1.Con,1,function(x){length(which(x))}) / BootstrapTimes
PVal_1.Vec = cbind(Con.Vec,PVal_1.sim)
PVal_2.sim = apply(Bootstrap_2.Con,1,function(x){length(which(x))}) / BootstrapTimes
PVal_2.Vec = cbind(Con.Vec,PVal_2.sim)
proc.time() - ptm





# ------------------------------------------------------------------------------------
# False Discovery Rate Resampling Adjustments 
# Yekutieli, D., & Benjamini, Y. (1999). Resampling-based false discovery rate controlling multiple test procedures for correlated test statistics. Journal of Statistical Planning and Inference, 82(1-2), 171-196. doi:10.1016/S0378-3758(99)00041-5
# Publishing, S. (2010). SAS/STAT 9. 22 User's Guide. The MIXED Procedure (Book Excerpt) (p. 228). SAS Institute.
# ------------------------------------------------------------------------------------

# different node space 
GetPValByFDR <- function(PVal.Vec){
	FDR_m = nrow(Con.Vec)
	Node.Space = unique(sort(Con.Vec[,2]-Con.Vec[,1]))

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
PVal_1.FDR.Vec = cbind(Con.Vec,GetPValByFDR(PVal_1.Vec))
PVal_2.FDR.Vec = cbind(Con.Vec,GetPValByFDR(PVal_2.Vec))

# ------------------------------------------------------------------------------------
CircosPickStr <- function(PValue.Vec,PickInd){
    if(length(PickInd)==0){return()}
    iCount = 1
    LinkStr = ""
    for (k in 1:length(PickInd)){
        i = PickInd[k]
        iRow = PValue.Vec[i,1]
        iCol = PValue.Vec[i,2]
        PositionStr1.1 = GetExpPointPosFromEquPointInd(iRow)
        PositionStr1.2 = GetExpPointPosFromEquPointInd(iRow)
        PositionStr2.1 = GetExpPointPosFromEquPointInd(iCol)
        PositionStr2.2 = GetExpPointPosFromEquPointInd(iCol)
        Field1 = paste("segdup",sprintf("%04d hs16", iCount),sep="")
        params = sprintf("id=%g",2,sep="")
        LinkStr = paste(LinkStr, paste(Field1," ",PositionStr2.1," ",PositionStr2.2," ",params,"\n",sep=""),sep="")
        LinkStr = paste(LinkStr, paste(Field1," ",PositionStr1.1," ",PositionStr1.2," ",params,"\n",sep=""),sep="")
        iCount = iCount + 1
    }
    return(LinkStr)
}


PS2Dist <- function(PS.Vec,PickInd,CutColInd){
	if(length(PickInd)==0){return()}
	
#	Sigma    = CutOff/3
	MaxPS_ij = max(PS.Vec[PickInd,3])
  MinPS_ij = min(PS.Vec[PS.Vec[PickInd,3]!=0,3])
  Mu      = BallDiameter 
  Sigma = (CutOff - Mu)/sqrt( 2*log(MaxPS_ij / MinPS_ij) )
	D_ij     = Mu + Sigma*sqrt(2*log(MaxPS_ij/PS.Vec[PickInd,3]))
	D_ij = sprintf("%.3f",D_ij)
	return(cbind(Con.Vec[PickInd,],D_ij))
}

Alpha = 0.05
AlphaStr = round(Alpha*100)
PickInd_1  = which(PVal_1.FDR.Vec[,3]<Alpha)
LinkStr_1  = CircosPickStr(PVal_1.FDR.Vec,PickInd_1)
Dist_1.Vec = PS2Dist(PS_1.Vec,PickInd_1,CutColInd)

PickInd_2  = which(PVal_2.FDR.Vec[,3]<Alpha)
LinkStr_2  = CircosPickStr(PVal_2.FDR.Vec,PickInd_2)
Dist_2.Vec = PS2Dist(PS_2.Vec,PickInd_2,CutColInd)

of_name = paste(DestDir,"GM12878.",SampleSize,".c",CutOff,".a",AlphaStr,".exp.pval.pick.txt",sep="")
write.table(file=of_name,
     LinkStr_1,quote=FALSE,row.names=FALSE,col.names=FALSE)
of_name = paste(DestDir,"GM12878.",SampleSize,".c",CutOff,".a",AlphaStr,".sis.node.dst",sep="")
write.table(file=of_name,
     Dist_1.Vec,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
of_name = paste(DestDir,"K562.",SampleSize,".c",CutOff,".a",AlphaStr,".exp.pval.pick.txt",sep="")
write.table(file=of_name,
     LinkStr_2,quote=FALSE,row.names=FALSE,col.names=FALSE)
of_name = paste(DestDir,"K562.",SampleSize,".c",CutOff,".a",AlphaStr,".sis.node.dst",sep="")
write.table(file=of_name,
     Dist_2.Vec,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")


# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------

