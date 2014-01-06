#!/usr/bin/R

rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)

if (Sys.info()["sysname"] == "Darwin"){
	CELL = "K562"
	# CELL = "GM12878"
	SampleSize = "p600_d6000"
	alpha =0.001
	Cutoff = 600
	LengthDiscard = 40000
}else{
	CELL = args[1]
	SampleSize = args[2]
	alpha = as.numeric(args[3])
	Cutoff = as.numeric(args[4])
	LengthDiscard = as.numeric(args[5])
}



PersistenceLength= as.numeric(sub("p(.*)_d.*","\\1",SampleSize))
if (PersistenceLength == 600){
    SegmentBPLength = 2100
}else if (PersistenceLength == 900){
    SegmentBPLength = 3200
}else if (PersistenceLength == 1200){
    SegmentBPLength = 4300
}
#------------------------------------------------------------------
# Parameter
#------------------------------------------------------------------
SourcePath = "data/ENm008/"
# CELL = "GM12878"
ReferenceCountFile = paste("data/analysis/ENm008/beadstring.",SampleSize,".cut",Cutoff,".comb.txt",sep="")

# DestDir = "/tmp/tmp/"
DestDir = "tmp/"
#------------------------------------------------------------------
# Get MaxRead between GM12878 and K562
#------------------------------------------------------------------
ContactFrequencyFile = paste(SourcePath,"Nature_2010_5CFrequencyCountsMatrix_ENm008_GM12878.txt",
    sep="")
ContactFrequencyData_1 = read.table(ContactFrequencyFile)

ContactFrequencyFile = paste(SourcePath,"Nature_2010_5CFrequencyCountsMatrix_ENm008_K562.txt",
    sep="")
ContactFrequencyData_2 = read.table(ContactFrequencyFile)

#------------------------------------------------------------------
# Read Contact Frequency
#------------------------------------------------------------------
ContactFrequencyFile = paste(SourcePath,"Nature_2010_5CFrequencyCountsMatrix_ENm008_",CELL,".txt",
    sep="")
ContactFrequencyData = read.table(ContactFrequencyFile)

Row_Ind_Start = as.numeric(sub(".*chr16:(.*)-(.*)","\\1",rownames(ContactFrequencyData)))
Row_Ind_End = as.numeric(sub(".*chr16:(.*)-(.*)","\\2",rownames(ContactFrequencyData)))
Row_Ind_Diff = Row_Ind_Start[2:(length(Row_Ind_Start))] - Row_Ind_End[1:(length(Row_Ind_End)-1)] 

Col_Ind_Start = as.numeric(sub(".*chr16\\.(.*)\\.(.*)","\\1",colnames(ContactFrequencyData)))
Col_Ind_End = as.numeric(sub(".*chr16\\.(.*)\\.(.*)","\\2",colnames(ContactFrequencyData)))
Col_Ind_Diff = Col_Ind_Start[2:(length(Col_Ind_Start))] - Col_Ind_End[1:(length(Col_Ind_End)-1)] 



#------------------------------------------------------------------
# Reconstruct the reading count matrix, gap is considered.
#------------------------------------------------------------------
# Start from Col_Col_Start, Col_Ind_Start = 1
MinLength = min(Row_Ind_End-Row_Ind_Start,Col_Ind_End - Col_Ind_Start)

Row_Vec = cbind(Row_Ind_Start,Row_Ind_End)
Col_Vec = cbind(Col_Ind_Start,Col_Ind_End)
All_Vec = rbind(Row_Vec,Col_Vec)
All_Vec_Sort = sort(All_Vec[,1],index.return=T)
All_Vec = All_Vec[All_Vec_Sort$ix,]
All_Vec_Diff = (All_Vec[,1])[2:nrow(All_Vec)] - (All_Vec[,2])[1:(nrow(All_Vec)-1)]

# expand the gap sequence length less than min seq length
MinInd = which(All_Vec_Diff <MinLength & All_Vec_Diff != 0)
All_Vec[MinInd,2] = All_Vec[MinInd,2] + All_Vec_Diff[MinInd]
All_Vec_Diff[MinInd] = 0

# insert the gap and resort, get segment interval
Gap_Vec = cbind((All_Vec[,2])[1:(nrow(All_Vec)-1)], (All_Vec[,2])[1:(nrow(All_Vec)-1)] + All_Vec_Diff)
All_Vec = rbind(All_Vec,Gap_Vec[which(All_Vec_Diff!=0),])
All_Vec_Sort = sort(All_Vec[,1],index.return=T)
All_Vec = All_Vec[All_Vec_Sort$ix,]

nSize = nrow(All_Vec)
IntervalPoints = unique(sort(c(All_Vec)))

#------------------------------
# Output Segemntindex NodeStart and NodeEnd
#------------------------------
# of_name = paste(DestDir,CELL,".SegInd.Node_StartEnd.txt",sep="")
# write.table(file=of_name,
#     cbind(1:nrow(All_Vec),All_Vec),quote=FALSE,row.names=FALSE,col.names=FALSE)

#----------------------------------------------------------------------------
# GM12878 vector (i, j , count)
GetCFVec = function(ContactFrequencyData, LengthDiscard=0){
	CF.Vec.tmp = NULL
	for ( i in 1:nrow(ContactFrequencyData)){
	    for (j in 1:ncol(ContactFrequencyData)){
	        r1 = findInterval(Row_Vec[i,1], IntervalPoints)
	        r2 = findInterval(Row_Vec[i,2]-1, IntervalPoints)
	        c1 = findInterval(Col_Vec[j,1], IntervalPoints)
	        c2 = findInterval(Col_Vec[j,2]-1, IntervalPoints)
	        if (r1 == r2 & c1 == c2 ){
							if(r1 < c1){
								ijCon = c(r1,c1, ContactFrequencyData[i,j])
							}else{
								ijCon = c(c1,r1, ContactFrequencyData[i,j])
							}
	            CF.Vec.tmp = rbind(CF.Vec.tmp,ijCon)
	        }
	    }
	}
	row.names(CF.Vec.tmp) = NULL
	# remove neighbors
	CF.Vec.tmp = CF.Vec.tmp[which(CF.Vec.tmp[,2]-CF.Vec.tmp[,1]!=1),]
	# remove length less than x kbp, segment head position - segment tail position
	if (LengthDiscard != 0){
		CF.Vec.tmp = CF.Vec.tmp[which(IntervalPoints[CF.Vec.tmp[,2]] - IntervalPoints[CF.Vec.tmp[,1]+1] >LengthDiscard),]
	}
	
	# sort first and second
	CF.Vec.tmp = CF.Vec.tmp[sort.int(CF.Vec.tmp[,1],index.return=T)$ix,]
	CF.Vec.tmp = CF.Vec.tmp[sort.int(CF.Vec.tmp[,2],index.return=T)$ix,]
	return(CF.Vec.tmp)
}

CF_1.Vec = GetCFVec(ContactFrequencyData_1,LengthDiscard)
# plot(CF_1.Vec[,3])
CF_2.Vec = GetCFVec(ContactFrequencyData_2,LengthDiscard)
# plot(CF_2.Vec[,3])

# plot(CF_2.Vec[which(IntervalPoints[CF_2.Vec[,2]] - IntervalPoints[CF_2.Vec[,1]] > 10000),3])
#----------------------------------------------------------------------------
# circos output histogram
# MPG is housekeeping gene
HS48.1.SegInd = findInterval(91497,IntervalPoints)
HS48.2.SegInd = findInterval(95256,IntervalPoints)
HS40SegInd = findInterval(100530,IntervalPoints)
HS10SegInd = findInterval(131220,IntervalPoints)
A2SegInd = findInterval(147782,IntervalPoints)
MPGSegInd = findInterval(65723,IntervalPoints)
POLSegInd = findInterval(29779,IntervalPoints)

#----------------------------------------------------------------------------
# current cell matrix
if (CELL == "GM12878"){
    CF.Vec = CF_1.Vec
}else{ # CELL is K562
    CF.Vec = CF_2.Vec
}
# sort first and second
CF.Vec = CF.Vec[sort.int(CF.Vec[,1],index.return=T)$ix,]
CF.Vec = CF.Vec[sort.int(CF.Vec[,2],index.return=T)$ix,]

# record index arr fo contact frequence do not equal to 0
ConInd.Vec.Ind = which(CF.Vec[,3]!=0,arr.ind=T)
ConInd.Vec = CF.Vec[ConInd.Vec.Ind,1:2]


CF.Vec = CF.Vec[ConInd.Vec.Ind,]
CF.SumCount = sum(CF.Vec[,3])
CF.Vec.Scale = CF.Vec
CF.Vec.Scale[,3] = CF.Vec[,3]/CF.SumCount


#------------------------------
# Persistence length
#------------------------------
# Dekker, J. (2008). Mapping in vivo chromatin interactions in yeast suggests an extended chromatin fiber with regional variation in compaction The Journal of biological chemistry, 283(50), 34532-34540. doi:10.1074/jbc.M806479200
# the mass density 28 nm/kb
# the persistence length is 58-118 nm, the mean is 88 nm
#------------------------------
# Fussner, E., Ching, R. W., & Bazett-Jones, D. P. (2011). Living without 30nm chromatin fibers. Trends in biochemical sciences, 36(1), 1-6. doi:10.1016/j.tibs.2010.09.002
# 10 nm chromatin fiber
#------------------------------
# 90 nm / 28 nm/kb = 3214 bp \approx 3200 bp
# 10 nm / 28 nm/kb = 357  bp \approx 350  bp
# One segment length 88 nm - 10 nm = 78 nm, containing 3000-350 bp
# the diameter of ball is 10 nm. the ball contains 350 bp
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


#------------------------------------------------------------------
# Read reference data, and merge to experimental segment
#------------------------------------------------------------------
Read.RF.Data = function(SegStartEnd.Forward.List,RFName,LengthDiscard=0){
    RFData = read.table(RFName)
    rowSize = length(SegStartEnd.Forward.List)
    segSize = max(SegStartEnd.Forward.List[[rowSize]])
    RF.Vec = matrix(0,nrow=(rowSize*(rowSize-1)/2),ncol=3)

    RF.Row.Ind.List = list()
    RF.Col.Ind.List = list()
    for (i in 1:segSize){
        MatchInd = which(RFData[,1]==i)
        if(length(MatchInd)>0){
            RF.Row.Ind.List[[i]] = MatchInd
        }
        else{
            RF.Row.Ind.List[[i]] = 0
        }
        
        MatchInd = which(RFData[,2]==i)
        if(length(MatchInd)>0){
            RF.Col.Ind.List[[i]] = MatchInd
        }
        else{
            RF.Col.Ind.List[[i]] = 0
        }
    }


    iCount = 1
    for (i in 1:(rowSize-1)){
        RFData.Seg.Row.Vec = SegStartEnd.Forward.List[[i]]
        RFData.RowInd = c()
        for (ii in 1:length(RFData.Seg.Row.Vec)){
            RFData.RowInd = c(RFData.RowInd,RF.Row.Ind.List[[RFData.Seg.Row.Vec[ii]]])
        }
        RFData.RowInd = unique(sort(RFData.RowInd))

        for (j in (i+1):rowSize){
            RFData.Seg.Col.Vec = SegStartEnd.Forward.List[[j]]
            RFData.ColInd = c()
            for (jj in 1:length(RFData.Seg.Col.Vec)){
                RFData.ColInd = c(RFData.ColInd,RF.Col.Ind.List[[RFData.Seg.Col.Vec[jj]]])
            }
            RFData.ColInd = unique(sort(RFData.ColInd))
            
            RFData.Ind = intersect(RFData.RowInd, RFData.ColInd)

            if (length(RFData.Ind)>0){
                RF.Vec[iCount,] = c(i,j,max(RFData[RFData.Ind,3]))
            }else{
                RF.Vec[iCount,] = c(i,j,0)
            }
            iCount = iCount +1
        }
    }

		# remove neighbor contacts
		RF.Vec = RF.Vec[RF.Vec[,2]-RF.Vec[,1]!=1,]
		# remove length less than x kbp, segment head position - segment tail position
		if (LengthDiscard != 0){
			RF.Vec = RF.Vec[which(IntervalPoints[RF.Vec[,2]] - IntervalPoints[RF.Vec[,1]+1] >LengthDiscard),]
		}

		# sort first and second
		RF.Vec = RF.Vec[sort.int(RF.Vec[,1],index.return=T)$ix,]
		RF.Vec = RF.Vec[sort.int(RF.Vec[,2],index.return=T)$ix,]

    return (RF.Vec)
}

RF.Vec = Read.RF.Data(SegStartEnd.Forward.List,ReferenceCountFile,LengthDiscard)
# #------------------------------------------------------------------
# #------------------------------------------------------------------

RF.Vec.Scale = RF.Vec
RF.SumCount = sum(RF.Vec[,3])
RF.Vec.Scale[,3] = RF.Vec.Scale[,3]/RF.SumCount

#------------------------------------------------------------------
# Ensemble contact for circos
#------------------------------------------------------------------

#------------------------------------------------------------------
# Propensity
#------------------------------------------------------------------
PS.Vec = cbind(ConInd.Vec,rep(0,nrow(ConInd.Vec)))

for ( i in 1:nrow(ConInd.Vec)){
    iRow = ConInd.Vec[i,1]
    iCol = ConInd.Vec[i,2]
    RFVecInd = which(RF.Vec[,1]==iRow & RF.Vec[,2]==iCol)
    
    if (length(RFVecInd) == 0 | RF.Vec[RFVecInd,3] == 0){
        PS.Vec[i,3] = CF.Vec.Scale[i,3] / (min(RF.Vec[RF.Vec[,3]!=0,3])/RF.SumCount)
    }else{
        PS.Vec[i,3] = CF.Vec.Scale[i,3] / RF.Vec.Scale[RFVecInd,3]
    }
}
PS.Max = max(PS.Vec[,3])


#------------------------------------------------------------------

#------------------------------------------------------------------
# different alpha
#------------------------------------------------------------------

PSToSeg = function(PS.Vec, PS.Vec.Sort.Ind, PickSize){
	if (PickSize==0){return(NULL)}
	#--------------------
	PSVal.Output = PS.Vec[PS.Vec.Sort.Ind[1:PickSize],,drop=FALSE]
	# sort first segment
	PSVal.Output = PSVal.Output[sort.int(PSVal.Output[,1],index.return=T)$ix,,drop=FALSE]
	# sort second segment
	PSVal.Output = PSVal.Output[sort.int(PSVal.Output[,2],index.return=T)$ix,,drop=FALSE]

	Seg.PVal = c()
	for (i in 1:nrow(PSVal.Output)){
		iRowInd.Exp = PSVal.Output[i,1]
		iColInd.Exp = PSVal.Output[i,2]
		PS.Exp = PSVal.Output[i,3]

		iRowInd.Seg = SegStartEnd.Forward.List[[iRowInd.Exp]]
		iColInd.Seg = SegStartEnd.Forward.List[[iColInd.Exp]]
		SegIJ = expand.grid(iRowInd.Seg, iColInd.Seg)
		kInd = apply(SegIJ,1,function(x){x[1]!=x[2]})
		if (length(which(kInd))==0){next}
		SegIJ = SegIJ[ kInd ,]
		SegIJ = cbind(SegIJ, PS.Exp)
		Seg.PVal = rbind(Seg.PVal, SegIJ)
	}

	Seg.PVal = Seg.PVal[sort.int(Seg.PVal[,1],index.return=T)$ix,]
	Seg.PVal = Seg.PVal[sort.int(Seg.PVal[,2],index.return=T)$ix,]
	
	# calculate distance
	Sigma = Cutoff/3
	MaxP_ij = max(Seg.PVal[,3])
	Mu = 100+2.4
	D_ij = Mu + Sigma*sqrt(2*log(MaxP_ij/Seg.PVal[,3]))
	Seg.PVal[,4] = sprintf("%.3f",D_ij)
	Seg.PVal[,3] = sprintf("%.3f",Seg.PVal[,3])

	rownames(Seg.PVal) = 1:nrow(Seg.PVal)
	return(Seg.PVal)
}

CircosPickStr = function(PValue.Vec,PValue.Vec.Sort.Ind,PickSize,Type,iCount){
    if(PickSize==0){return()}
    # iCount = 1
    LinkStr = ""
    for (k in 1:PickSize){
        i = PValue.Vec.Sort.Ind[k]
        iRow = PValue.Vec[i,1]
        iCol = PValue.Vec[i,2]
        PositionStr1.1 = IntervalPoints[iRow]
        PositionStr1.2 = IntervalPoints[iRow+1]
        PositionStr2.1 = IntervalPoints[iCol]
        PositionStr2.2 = IntervalPoints[iCol+1]
        Field1 = paste("segdup",sprintf("%04d hs16", iCount),sep="")
        params = sprintf("id=%g",Type,sep="")
        LinkStr = paste(LinkStr, paste(Field1," ",PositionStr2.1," ",PositionStr2.2," ",params,sep=""),sep="\n")
        LinkStr = paste(LinkStr, paste(Field1," ",PositionStr1.1," ",PositionStr1.2," ",params,sep=""),sep="\n")
        iCount = iCount + 1
    }
    return(LinkStr)
}


UniqSeg = function(Seg.PVal){
    if(is.null(Seg.PVal)){return(NULL)}
    DupInd = which(duplicated(Seg.PVal[,1:2]))
    setdiff(1:nrow(Seg.PVal),c(DupInd-1,DupInd))
    PrvAftInd = unique(sort(c(DupInd-1,DupInd)))

    Uniq1.PVal = Seg.PVal[setdiff(1:nrow(Seg.PVal),PrvAftInd),]

    Uniq2.PVal = c()
    for (i in 1:length(DupInd)){
        SegDupInd = DupInd[i]
        AllDupInd = which(Seg.PVal[,1] == Seg.PVal[SegDupInd,1] & Seg.PVal[,2] == Seg.PVal[SegDupInd,2] )
        if (length(unique(Seg.PVal[AllDupInd,3])) == 1 & 
            (length(which(Uniq2.PVal[,1]==Seg.PVal[SegDupInd,1] & Uniq2.PVal[,2]==Seg.PVal[SegDupInd,2])))==0
            ){
            Uniq2.PVal = rbind(Uniq2.PVal,Seg.PVal[SegDupInd,])
        }
    }

    Uniq.PVal = rbind(Uniq1.PVal,Uniq2.PVal)
    Uniq.PVal = Uniq.PVal[sort.int(Uniq.PVal[,1],index.return=T)$ix,]
    Uniq.PVal = Uniq.PVal[sort.int(Uniq.PVal[,2],index.return=T)$ix,]
    return(Uniq.PVal)
}
# ------------------------------------------------------------------

PickFile=paste("tmp/",CELL,".",SampleSize,".cut", Cutoff, ".a",alpha,".l",LengthDiscard,".exp.pval.pick.txt" ,sep="")
Pick.Raw = read.table(PickFile)
Pick.Size = nrow(Pick.Raw)/2
Pick.Ind.Vec = matrix(0,nrow=Pick.Size,ncol=2)
for (i in 1:Pick.Size){
	Pick.Ind.Vec[i,] = c(which(IntervalPoints==Pick.Raw[2*(i-1)+2,3]),which(IntervalPoints==Pick.Raw[2*(i-1)+1,3]) )
}
Pick.Ind.Vec = cbind(Pick.Ind.Vec,apply(Pick.Ind.Vec,1,function(x){PS.Vec[PS.Vec[,1]==x[1] & PS.Vec[,2]==x[2],3]}))

of_name = paste("tmp/img/",CELL,".",SampleSize,".cut", Cutoff, ".a",alpha,".l",LengthDiscard,".hist.pdf" ,sep="")
pdf(of_name)
hist(Pick.Ind.Vec[,3],breaks=100,xlab="Propensity",main=sub(".hist.png","",basename(of_name)))
dev.off()

#----------------------------------------
# circos of Propensity
#----------------------------------------

CircosPick_PSStr = function(Pick.Ind.Vec_,iCount){
		PickSize = nrow(Pick.Ind.Vec_)
    if(PickSize==0){return()}
		MaxPS = max(Pick.Ind.Vec_[,3])
    # iCount = 1
    LinkStr = ""
    for (k in 1:PickSize){
        iRow = Pick.Ind.Vec_[k,1]
        iCol = Pick.Ind.Vec_[k,2]
        PositionStr1.1 = IntervalPoints[iRow]
        PositionStr1.2 = IntervalPoints[iRow+1]
        PositionStr2.1 = IntervalPoints[iCol]
        PositionStr2.2 = IntervalPoints[iCol+1]
        Field1 = paste("segdup",sprintf("%04d hs16", iCount),sep="")
				freq = sprintf("thickness=%.3f",Pick.Ind.Vec_[k,3])
				maxfreq = sprintf("id=%.2f",MaxPS,sep="")
				params = paste(freq,",",maxfreq,sep="")
        LinkStr = paste(LinkStr, paste(Field1," ",PositionStr2.1," ",PositionStr2.2," ",params,sep=""),sep="\n")
        LinkStr = paste(LinkStr, paste(Field1," ",PositionStr1.1," ",PositionStr1.2," ",params,sep=""),sep="\n")
        iCount = iCount + 1
    }
    return(LinkStr)
}

PS_threshold_List = c(5,10,20,40,60)

for (i in 1:length(PS_threshold_List)){
	PS_threshold = PS_threshold_List[i]

	iCount = 1
	Pick.Cut.Ind.Vec = Pick.Ind.Vec[Pick.Ind.Vec[,3]>PS_threshold,,drop=FALSE]
	cat("PickCutSize= ", nrow(Pick.Cut.Ind.Vec), "\n")
	if(nrow(Pick.Cut.Ind.Vec)==0){
		break
	}
	LinkStr = CircosPick_PSStr(Pick.Cut.Ind.Vec,iCount)
	of_circos = paste("tmp/img/",CELL,".",SampleSize,".cut", Cutoff, ".a",alpha,".l",LengthDiscard,".p",PS_threshold,".prop.txt" ,sep="")
	write.table(file=of_circos,LinkStr,quote=FALSE,row.names=FALSE,col.names=FALSE)
}



# ------------------------------------------------------------------
# ------------------------------------------------------------------
# ------------------------------------------------------------------
# ------------------------------------------------------------------
# ------------------------------------------------------------------
# ------------------------------------------------------------------

