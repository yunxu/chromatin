#!/usr/bin/R
###################################################################
# Get contact index according to segment (starting and end)
###################################################################

rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
fileid = args[1]
#------------------------------------------------------------------
# Parameter
#------------------------------------------------------------------
SourcePath = "data/ENm008/"
CELL = "GM12878"
#CELL = "K562"
ReferenceCountFile = paste("data/analysis/ENm008/64.comb.txt",sep="")
BinSize = 60
DestDir = "/dump/tmp/"

#------------------------------------------------------------------
# Get MaxRead between GM12878 and K562
#------------------------------------------------------------------
ContactFrequencyFile = paste(SourcePath,"Nature_2010_5CFrequencyCountsMatrix_ENm008_GM12878.txt",
	sep="")
ContactFrequencyData_1 = read.table(ContactFrequencyFile)
Ratio_1 = max(ContactFrequencyData_1)/sum(ContactFrequencyData_1)
# > Ratio_1
# [1] 0.03182158
# > sum(ContactFrequencyData_1)
# [1] 182989
# > max(ContactFrequencyData_1)
# [1] 5823

ContactFrequencyFile = paste(SourcePath,"Nature_2010_5CFrequencyCountsMatrix_ENm008_K562.txt",
	sep="")
ContactFrequencyData_2 = read.table(ContactFrequencyFile)
Ratio_2 = max(ContactFrequencyData_2)/sum(ContactFrequencyData_2)
# > Ratio_2
# [1] 0.1037235
# > sum(ContactFrequencyData_2)
# [1] 131947
# > max(ContactFrequencyData_2)
# [1] 13686
# MaxCount = max(ContactFrequencyData_1,ContactFrequencyData_2)
Ratio_Max = max(Ratio_1,Ratio_2)
#------------------------------------------------------------------
# Read Contact Frequency
#------------------------------------------------------------------

ContactFrequencyFile = paste(SourcePath,"Nature_2010_5CFrequencyCountsMatrix_ENm008_",CELL,".txt",
	sep="")
ContactFrequencyData = read.table(ContactFrequencyFile)
MaxCount = sum(ContactFrequencyData)
# MaxCount = max(ContactFrequencyData_1,ContactFrequencyData_2)



Row_Ind_Start = as.numeric(sub(".*chr16:(.*)-(.*)","\\1",rownames(ContactFrequencyData)))
Row_Ind_End = as.numeric(sub(".*chr16:(.*)-(.*)","\\2",rownames(ContactFrequencyData)))
Row_Ind_Diff = Row_Ind_Start[2:(length(Row_Ind_Start))] - Row_Ind_End[1:(length(Row_Ind_End)-1)] 

Col_Ind_Start = as.numeric(sub(".*chr16\\.(.*)\\.(.*)","\\1",colnames(ContactFrequencyData)))
Col_Ind_End = as.numeric(sub(".*chr16\\.(.*)\\.(.*)","\\2",colnames(ContactFrequencyData)))
Col_Ind_Diff = Col_Ind_Start[2:(length(Col_Ind_Start))] - Col_Ind_End[1:(length(Col_Ind_End)-1)] 


#------------------------------------------------------------------
# Output for circos highlight
#------------------------------------------------------------------
# Row_Str = paste("hs16",Row_Ind_Start,Row_Ind_End,rep(c("fill_color=lgreen","fill_color=green")))
# Col_Str = paste("hs16",Col_Ind_Start,Col_Ind_End,rep(c("fill_color=lblue","fill_color=blue")))
# 
# of_circos = paste(DestDir,"raw.segment.for.highlight.txt",sep="")
# write.table(file=of_circos,Row_Str,quote=FALSE,row.names=FALSE,col.names=FALSE)
# of_circos = paste(DestDir,"raw.segment.rev.highlight.txt",sep="")
# write.table(file=of_circos,Col_Str,quote=FALSE,row.names=FALSE,col.names=FALSE)
# 
# Row_Sep_Str = paste(
# 	"hs16", Row_Ind_End[which(Row_Ind_Diff==0)]-2, Row_Ind_End[which(Row_Ind_Diff==0)]+2)
# of_circos = paste(DestDir,"raw.segment.for.sep.highlight.txt",sep="")
# write.table(file=of_circos,Row_Sep_Str,quote=FALSE,row.names=FALSE,col.names=FALSE)
# 
# Col_Sep_Str = paste(
# 	"hs16", Col_Ind_End[which(Col_Ind_Diff==0)]-2, Col_Ind_End[which(Col_Ind_Diff==0)]+2)
# of_circos = paste(DestDir,"raw.segment.rev.sep.highlight.txt",sep="")
# write.table(file=of_circos,Col_Sep_Str,quote=FALSE,row.names=FALSE,col.names=FALSE)


# iCount = 1
# LinkStr = ""
# for (i in 1:(length(Row_Ind_Start))){
# 	for (j in 1:(length(Col_Ind_Start))){
# 		if(ContactFrequencyData[i,j]==0){next}
# 		PositionStr1.1 = Row_Ind_Start[i]
# 		PositionStr1.2 = Row_Ind_End[i]
# 		PositionStr2.1 = Col_Ind_Start[j]
# 		PositionStr2.2 = Col_Ind_End[j]
# 		Field1 = paste("segdup",sprintf("%04d hs16", iCount),sep="")
# 		freq = sprintf("thickness=%.3f",ContactFrequencyData[i,j],sep="")
# 		# maxfreq = sprintf("id=%.3f",MaxCount * Ratio_Max,sep="")
# 		# maxfreq = sprintf("id=%.3f", max(ContactFrequencyData),sep="")
# 		maxfreq = sprintf("id=%.3f",MaxCount,sep="")
# 		params = paste(freq,",",maxfreq,sep="")
# 		LinkStr = paste(LinkStr, paste(Field1," ",PositionStr2.1," ",PositionStr2.2," ",params,sep=""),sep="\n")
# 		LinkStr = paste(LinkStr, paste(Field1," ",PositionStr1.1," ",PositionStr1.2," ",params,sep=""),sep="\n")
# 		iCount = iCount + 1
# 	}
# }
# 
# of_circos = paste(DestDir,"contactfreq.",CELL,".raw.txt",sep="")
# write.table(file=of_circos,LinkStr,quote=FALSE,row.names=FALSE,col.names=FALSE)

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

CF.Matrix = matrix(0,nrow=nSize,ncol=nSize)
for ( i in 1:nrow(ContactFrequencyData)){
	for (j in 1:ncol(ContactFrequencyData)){
		r1 = findInterval(Row_Vec[i,1], IntervalPoints)
		r2 = findInterval(Row_Vec[i,2]-1, IntervalPoints)
		c1 = findInterval(Col_Vec[j,1], IntervalPoints)
		c2 = findInterval(Col_Vec[j,2]-1, IntervalPoints)
		if (r1 == r2 & c1 == c2 ){
			CF.Matrix[r1,c1] = CF.Matrix[c1,r1] = ContactFrequencyData[i,j]
		}
	}
}

# recode index arr fo contact frequence do not equal to 0
ConInd_Vec = which(CF.Matrix!=0,arr.ind=T)
# only get upper right index
ConInd_Vec = ConInd_Vec[which(ConInd_Vec[,1]<ConInd_Vec[,2]),]
# scale matrix
CF.Matrix.Scale = CF.Matrix/MaxCount

#------------------------------------------------------------------
# Find interval starting index and end index of same sapce
#------------------------------------------------------------------
GetStartEnd = function(SegPoints){
	NumPoints = length(SegPoints)
	NumSeg = NumPoints - 1
	MinSpace = min(SegPoints[2:NumPoints] - SegPoints[1:(NumPoints-1)])
	EquSpacePoints = seq(1,max(SegPoints),by=MinSpace)
  PointsInd = findInterval(EquSpacePoints,SegPoints)
	StartEndInd.Vec = matrix(0,nrow=NumSeg,ncol=2)
	for (i in 1:NumSeg){
		StartEndInd.Vec[i,] = c(min(which(PointsInd==i)),max(which(PointsInd==i)))
	}

	ConStartEndInd = unique(sort(StartEndInd.Vec))
	ConStartEndSeg = cbind(EquSpacePoints, EquSpacePoints+MinSpace-1)
	return(list(ConStartEndInd=ConStartEndInd,ConStartEndSeg=ConStartEndSeg))
}
ConStartEnd = GetStartEnd(IntervalPoints)

#--------------------
# output contact index file
#--------------------
# of_name = paste(DestDir,"ENm008_Cont_Index.txt",sep="")
# write.table(file=of_name,ConStartEnd$ConStartEndInd,quote=FALSE,row.names=FALSE,col.names=FALSE)

#--------------------
# output same segment starting and end point
#--------------------
# of_name = paste(DestDir,"ENm008_StartEnd_SameSeg.txt",sep="")
# write.table(file=of_name,ConStartEnd$ConStartEndSeg,quote=FALSE,row.names=FALSE,col.names=FALSE)

#------------------------------------------------------------------
# Read reference data, and merge to experimental segment
#------------------------------------------------------------------
RFData = read.table(ReferenceCountFile)
Seg1 = findInterval((RFData[,1]-1)*BinSize+1,IntervalPoints)
Seg2 = findInterval((RFData[,2]-1)*BinSize+1,IntervalPoints)
RF.Seg.Data = cbind(Seg1,Seg2,RFData[,3])

RF.Vec = cbind(ConInd_Vec,rep(0,nrow(ConInd_Vec)))
for (i in 1:nrow(ConInd_Vec)){
	iRow = ConInd_Vec[i,1]
	iCol = ConInd_Vec[i,2]
	MatchInd = which(RF.Seg.Data[,1]==iRow & RF.Seg.Data[,2]==iCol)
	if(length(MatchInd)>0){
		RF.Vec[i,3] = max(RF.Seg.Data[MatchInd,3])
	}
}
RF.MaxCount = max(RF.Seg.Data)
RF.Vec.Scale = RF.Vec
RF.Vec.Scale[,3] = RF.Vec.Scale[,3]/ RF.MaxCount

#------------------------------------------------------------------
# Ensemble contact for circos
#------------------------------------------------------------------

# iCount = 1
# LinkStr = ""
# for (i in 1:nrow(RF.Vec)){
# 	iRow = RF.Vec[i,1]
# 	iCol = RF.Vec[i,2]
# 	if(RF.Vec[i,3]==0){next}
# 	PositionStr1.1 = IntervalPoints[iRow]
# 	PositionStr1.2 = IntervalPoints[iRow+1]
# 	PositionStr2.1 = IntervalPoints[iCol]
# 	PositionStr2.2 = IntervalPoints[iCol+1]
# 	Field1 = paste("segdup",sprintf("%04d hs16", iCount),sep="")
# 	freq = sprintf("thickness=%.2f",RF.Vec[i,3])
# 	maxfreq = sprintf("id=%.2f",RF.MaxCount,sep="")
# 	params = paste(freq,",",maxfreq,sep="")
# 	LinkStr = paste(LinkStr, paste(Field1," ",PositionStr2.1," ",PositionStr2.2," ",params,sep=""),sep="\n")
# 	LinkStr = paste(LinkStr, paste(Field1," ",PositionStr1.1," ",PositionStr1.2," ",params,sep=""),sep="\n")
# 	iCount = iCount + 1
# }
# 
# of_circos = paste(DestDir,"contactfreq.ensemble.txt",sep="")
# write.table(file=of_circos,LinkStr,quote=FALSE,row.names=FALSE,col.names=FALSE)

#------------------------------------------------------------------
# Propensity
#------------------------------------------------------------------
PS.Vec = cbind(ConInd_Vec,rep(0,nrow(ConInd_Vec)))

for ( i in 1:nrow(ConInd_Vec)){
	iRow = ConInd_Vec[i,1]
	iCol = ConInd_Vec[i,2]
	if (RF.Vec[i,3] == 0){
		PS.Vec[i,3] = CF.Matrix.Scale[iRow,iCol] / (1/RF.MaxCount)
	}else{
		PS.Vec[i,3] = CF.Matrix.Scale[iRow,iCol] / RF.Vec.Scale[i,3]
	}
}
PS.Max = max(PS.Vec[,3])
# PS.Max = 322.68
# PS.Max = 4314.39

# iCount = 1
# LinkStr = ""
# for (i in 1:nrow(ConInd_Vec)){
# 	iRow = ConInd_Vec[i,1]
# 	iCol = ConInd_Vec[i,2]
# 	if(PS.Vec[i,3]==0){next}
# 	PositionStr1.1 = IntervalPoints[iRow]
# 	PositionStr1.2 = IntervalPoints[iRow+1]
# 	PositionStr2.1 = IntervalPoints[iCol]
# 	PositionStr2.2 = IntervalPoints[iCol+1]
# 	Field1 = paste("segdup",sprintf("%04d hs16", iCount),sep="")
# 	freq = sprintf("thickness=%.2f",PS.Vec[i,3])
# 	maxfreq = sprintf("id=%.2f",PS.Max,sep="")
# 	params = paste(freq,",",maxfreq,sep="")
# 	LinkStr = paste(LinkStr, paste(Field1," ",PositionStr2.1," ",PositionStr2.2," ",params,sep=""),sep="\n")
# 	LinkStr = paste(LinkStr, paste(Field1," ",PositionStr1.1," ",PositionStr1.2," ",params,sep=""),sep="\n")
# 	iCount = iCount + 1
# }
# 
# of_circos = paste(DestDir,"prop.",CELL,".txt",sep="")
# write.table(file=of_circos,LinkStr,quote=FALSE,row.names=FALSE,col.names=FALSE)


#------------------------------------------------------------------
# Permutation test for P-value
#------------------------------------------------------------------

# PermutaionTest = function(xcount,xtotal,ycount,ytotal){
# 	m = 1000
# 	nx = xtotal
# 	ny = ytotal
# 	nz = nx+ny
# 	zcount = xcount + ycount
# 	MeanDiff = xcount/nx - ycount/ny
# 	DiffVec = rep(0,m)
# 	for( i in 1:m){
# 		sampleind = sample(nz,zcount)
# 		xcount.rnd = length(which(sampleind < nx))
# 		ycount.rnd = zcount-xcount.rnd
# 		DiffVec[i] = xcount.rnd/nx - ycount.rnd/ny
# 	}
# 	#p-value for two-side test
# 	pvalue = sum(abs(DiffVec) > abs(MeanDiff)) /m
# 	# if (pvalue ==0){ pvalue=1/m }
# }
# 
# 
# Sign.Vec = cbind(ConInd_Vec,rep(0,nrow(ConInd_Vec)))
# for (i in 1:nrow(ConInd_Vec)){
# 	Diff = CF.Matrix[ConInd_Vec[i,1],ConInd_Vec[i,2]]/MaxCount - RF.Vec[i,3]/RF.MaxCount
# 	if (Diff > 0){
# 		Sign.Vec[i,3] = 1
# 	}else{
# 		Sign.Vec[i,3] = -1		
# 	}
# }
# 
# PValue.Vec = cbind(ConInd_Vec,rep(0,nrow(ConInd_Vec)))
# for (i in 1:nrow(ConInd_Vec)){
# 	cat(i,"/",nrow(ConInd_Vec),"\n")
# 	PValue.Vec[i,3] = PermutaionTest(CF.Matrix[ConInd_Vec[i,1],ConInd_Vec[i,2]],MaxCount, RF.Vec[i,3], RF.MaxCount)
# }
# 
# P.Cut = 0.05
# P.SigInd = p.adjust(PValue.Vec[,3],method="fdr") <= P.Cut
# 
# iCount = 1
# LinkStr = ""
# for (i in 1:nrow(PValue.Vec)){
# 	iRow = PValue.Vec[i,1]
# 	iCol = PValue.Vec[i,2]
# 	if(P.SigInd[i]==F | Sign.Vec[i,3]==-1){next}
# 	PositionStr1.1 = IntervalPoints[iRow]
# 	PositionStr1.2 = IntervalPoints[iRow+1]
# 	PositionStr2.1 = IntervalPoints[iCol]
# 	PositionStr2.2 = IntervalPoints[iCol+1]
# 	Field1 = paste("segdup",sprintf("%04d hs16", iCount),sep="")
# 	freq = sprintf("thickness=%.2f",0.1)
# 	maxfreq = sprintf("id=%d",1,sep="")
# 	params = paste(freq,",",maxfreq,sep="")
# 	LinkStr = paste(LinkStr, paste(Field1," ",PositionStr2.1," ",PositionStr2.2," ",params,sep=""),sep="\n")
# 	LinkStr = paste(LinkStr, paste(Field1," ",PositionStr1.1," ",PositionStr1.2," ",params,sep=""),sep="\n")
# 	iCount = iCount + 1
# }
# 
# of_circos = paste(DestDir,"pval.pos.",CELL,".txt",sep="")
# write.table(file=of_circos,LinkStr,quote=FALSE,row.names=FALSE,col.names=FALSE)
# 
# iCount = 1
# LinkStr = ""
# for (i in 1:nrow(PValue.Vec)){
# 	iRow = PValue.Vec[i,1]
# 	iCol = PValue.Vec[i,2]
# 	if(P.SigInd[i]==F | Sign.Vec[i,3]==1){next}
# 	PositionStr1.1 = IntervalPoints[iRow]
# 	PositionStr1.2 = IntervalPoints[iRow+1]
# 	PositionStr2.1 = IntervalPoints[iCol]
# 	PositionStr2.2 = IntervalPoints[iCol+1]
# 	Field1 = paste("segdup",sprintf("%04d hs16", iCount),sep="")
# 	freq = sprintf("thickness=%.2f",0.1)
# 	maxfreq = sprintf("id=%d",1,sep="")
# 	params = paste(freq,",",maxfreq,sep="")
# 	LinkStr = paste(LinkStr, paste(Field1," ",PositionStr2.1," ",PositionStr2.2," ",params,sep=""),sep="\n")
# 	LinkStr = paste(LinkStr, paste(Field1," ",PositionStr1.1," ",PositionStr1.2," ",params,sep=""),sep="\n")
# 	iCount = iCount + 1
# }
# 
# of_circos = paste(DestDir,"pval.neg.",CELL,".txt",sep="")
# write.table(file=of_circos,LinkStr,quote=FALSE,row.names=FALSE,col.names=FALSE)


#------------------------------------------------------------------
# Bootstrap for p-value
#------------------------------------------------------------------
# rebuild the distribution
BootstrapTest = function(xcount,xtotal,ycount,ytotal){
	# i=2
	# xcount = CF.Matrix[ConInd_Vec[i,1],ConInd_Vec[i,2]]
	# xtotal = MaxCount
	# ycount = RF.Vec[i,3]
	# ytotal = RF.MaxCount	
	
	x_mean = xcount/xtotal
	y_mean = ycount/ytotal
	sampleid = sample(ytotal,ycount)
	ReadCount = rep(0,ytotal)
	ReadCount[sampleid] = 1
	
	N_Trial = 1000
	Mean_Vec = rep(0,N_Trial)
	for(k in 1:N_Trial){
		# cat("\t",k,"/",N_Trial,"\n")
		BootStrapIDs =sample(ytotal,replace=T)
		Mean_Vec[k] = mean(ReadCount[BootStrapIDs])
	}
	Mean_Boot = mean(Mean_Vec)
	SE_Boot = sqrt((1/(N_Trial-1))* sum((Mean_Vec-Mean_Boot)^2))
	t_star = qt(c(0.025,0.975),df=N_Trial-1)
	CI_Boot = Mean_Boot+t_star*SE_Boot
	
	# less than or fall in or greater than in confidence interval
	l_e_g = findInterval(x_mean,CI_Boot)
	if (l_e_g == 2){
		P_Value = pt((x_mean-Mean_Boot)/SE_Boot,df=N_Trial-1, lower.tail=F)
	}else{
		P_Value = pt((x_mean-Mean_Boot)/SE_Boot,df=N_Trial-1, lower.tail=T)
	}
	return(list(l_e_g=l_e_g,P_Value=P_Value))
}

#--------------------
# Real P-Value
#--------------------

# PValue.Vec = cbind(ConInd_Vec,rep(0,nrow(ConInd_Vec)),rep(0,nrow(ConInd_Vec)))
# for (i in 1:nrow(ConInd_Vec)){
# 	cat(i,"/",nrow(ConInd_Vec),"\n")
#  	BootStrapTest = BootstrapTest(CF.Matrix[ConInd_Vec[i,1],ConInd_Vec[i,2]],MaxCount, RF.Vec[i,3], RF.MaxCount)
# 	PValue.Vec[i,3] = BootStrapTest$P_Value
# 	PValue.Vec[i,4] = BootStrapTest$l_e_g
# }
# 
# of_txt = paste(DestDir,"bootstrap_",CELL,"_pvalue.txt",sep="")
# write.table(file=of_txt,PValue.Vec,quote=FALSE,row.names=FALSE,col.names=FALSE)
# PValue.Vec = read.table(of_txt)
# 
# P.Cut = 0.05
# P.SigInd = p.adjust(PValue.Vec[,3],method="fdr") <= P.Cut

#--------------------
# FDR (False Discovery Rate)
#--------------------
N_FDR = 10
PValue.FDR.Matrix = matrix(0,nrow=nrow(ConInd_Vec),ncol=N_FDR)
for (i_fdr in 1:N_FDR){
	ShuffeledIDs = sample(nrow(ConInd_Vec))
	cat("fileid =",fileid,"i_fdr = ", i_fdr,"\n")
	for (i_con in 1:nrow(ConInd_Vec)){
    cat("\t i_con =",i_con,"/", "totalcon =",nrow(ConInd_Vec),"\n")
		BootStrapTest = BootstrapTest(CF.Matrix[ConInd_Vec[ShuffeledIDs[i_con],1],ConInd_Vec[ShuffeledIDs[i_con],2]],MaxCount, RF.Vec[i_con,3], RF.MaxCount)
		PValue.FDR.Matrix[i_con,i_fdr] = BootStrapTest$P_Value
	}
}

of_txt = paste(DestDir,"fdr_",CELL,"_",fileid,".txt",sep="")
write.table(file=of_txt,formatC(PValue.FDR.Matrix,format="e"),quote=F,col.names=F,row.names=F)
# iCount = 1
# LinkStr = ""
# for (i in 1:nrow(PValue.Vec)){
# 	iRow = PValue.Vec[i,1]
# 	iCol = PValue.Vec[i,2]
# 	if(P.SigInd[i]==F | PValue.Vec[i,4] != 2){next}
# 	PositionStr1.1 = IntervalPoints[iRow]
# 	PositionStr1.2 = IntervalPoints[iRow+1]
# 	PositionStr2.1 = IntervalPoints[iCol]
# 	PositionStr2.2 = IntervalPoints[iCol+1]
# 	Field1 = paste("segdup",sprintf("%04d hs16", iCount),sep="")
# 	freq = sprintf("thickness=%.2f",0.1)
# 	maxfreq = sprintf("id=%d",1,sep="")
# 	params = paste(freq,",",maxfreq,sep="")
# 	LinkStr = paste(LinkStr, paste(Field1," ",PositionStr2.1," ",PositionStr2.2," ",params,sep=""),sep="\n")
# 	LinkStr = paste(LinkStr, paste(Field1," ",PositionStr1.1," ",PositionStr1.2," ",params,sep=""),sep="\n")
# 	iCount = iCount + 1
# }
# 
# of_circos = paste(DestDir,"pval.pos.",CELL,".txt",sep="")
# write.table(file=of_circos,LinkStr,quote=FALSE,row.names=FALSE,col.names=FALSE)
# 
# iCount = 1
# LinkStr = ""
# for (i in 1:nrow(PValue.Vec)){
# 	iRow = PValue.Vec[i,1]
# 	iCol = PValue.Vec[i,2]
# 	if(P.SigInd[i]==F | PValue.Vec[i,4] != 0){next}
# 	PositionStr1.1 = IntervalPoints[iRow]
# 	PositionStr1.2 = IntervalPoints[iRow+1]
# 	PositionStr2.1 = IntervalPoints[iCol]
# 	PositionStr2.2 = IntervalPoints[iCol+1]
# 	Field1 = paste("segdup",sprintf("%04d hs16", iCount),sep="")
# 	freq = sprintf("thickness=%.2f",0.1)
# 	maxfreq = sprintf("id=%d",1,sep="")
# 	params = paste(freq,",",maxfreq,sep="")
# 	LinkStr = paste(LinkStr, paste(Field1," ",PositionStr2.1," ",PositionStr2.2," ",params,sep=""),sep="\n")
# 	LinkStr = paste(LinkStr, paste(Field1," ",PositionStr1.1," ",PositionStr1.2," ",params,sep=""),sep="\n")
# 	iCount = iCount + 1
# }
# 
# of_circos = paste(DestDir,"pval.neg.",CELL,".txt",sep="")
# write.table(file=of_circos,LinkStr,quote=FALSE,row.names=FALSE,col.names=FALSE)
