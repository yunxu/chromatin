#!/usr/bin/R
###################################################################
# Get contact index according to segment (starting and end)
###################################################################

rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
CELL = args[1]
SampleSize = args[2]
alpha = as.numeric(args[3])

# CELL = "K562"
# CELL = "GM12878"
# SampleSize = "640_r1100"
# SampleSize = "640_r6000"
# # SampleSize = "640_r60000"
# alpha =0.025
#------------------------------------------------------------------
# Parameter
#------------------------------------------------------------------
SourcePath = "data/ENm008/"
# CELL = "GM12878"
ReferenceCountFile = paste("data/analysis/ENm008/",SampleSize,".comb.txt",sep="")
BinSize = 60

DestDir = "/dump/temp/"
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

# 
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
# 
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
of_name = paste(DestDir,CELL,".SegInd.Node_StartEnd.txt",sep="")
write.table(file=of_name,
	cbind(1:nrow(All_Vec),All_Vec),quote=FALSE,row.names=FALSE,col.names=FALSE)



#----------------------------------------------------------------------------
# Get total neighbor count
GetNeighborCount = function(Matrix_){
	NeighborCount = 0
	for (i in 1:(nrow(Matrix_)-1)){
		NeighborCount = NeighborCount +Matrix_[i,i+1]
	}
	return (NeighborCount)
}

GetShortNeighborCount = function(Matrix_,IntervalPoints_){
	SegLen = IntervalPoints_[2:length(IntervalPoints_)] - IntervalPoints_[1:(length(IntervalPoints_)-1)]
	ShortInd = which(SegLen<5000)
	
	NeighborCount = 0
	for (i in 1:(length(ShortInd)-1)) {
		if (ShortInd[i+1]-ShortInd[i]==1){
			cat(sprintf("%d-%d %.2f\n", ShortInd[i], ShortInd[i+1], Matrix_[i,i+1]) )
		}
	}
}
#----------------------------------------------------------------------------
# GM12878 matrix
CF_1.Matrix = matrix(0,nrow=nSize,ncol=nSize)
for ( i in 1:nrow(ContactFrequencyData_1)){
	for (j in 1:ncol(ContactFrequencyData_1)){
		r1 = findInterval(Row_Vec[i,1], IntervalPoints)
		r2 = findInterval(Row_Vec[i,2]-1, IntervalPoints)
		c1 = findInterval(Col_Vec[j,1], IntervalPoints)
		c2 = findInterval(Col_Vec[j,2]-1, IntervalPoints)
		if (r1 == r2 & c1 == c2 ){
			CF_1.Matrix[r1,c1] = CF_1.Matrix[c1,r1] = ContactFrequencyData_1[i,j]
		}
	}
}

#----------------------------------------------------------------------------
# K562 matrix
CF_2.Matrix = matrix(0,nrow=nSize,ncol=nSize)
for ( i in 1:nrow(ContactFrequencyData_1)){
	for (j in 1:ncol(ContactFrequencyData_1)){
		r1 = findInterval(Row_Vec[i,1], IntervalPoints)
		r2 = findInterval(Row_Vec[i,2]-1, IntervalPoints)
		c1 = findInterval(Col_Vec[j,1], IntervalPoints)
		c2 = findInterval(Col_Vec[j,2]-1, IntervalPoints)
		if (r1 == r2 & c1 == c2 ){
			CF_2.Matrix[r1,c1] = CF_2.Matrix[c1,r1] = ContactFrequencyData_2[i,j]
		}
	}
}

#----------------------------------------------------------------------------
# GetShortNeighborCount(CF_1.Matrix,IntervalPoints)
# GetShortNeighborCount(CF_2.Matrix,IntervalPoints)


# circos output histogram
# MPG is housekeeping gene
HS40SegInd = findInterval(100530,IntervalPoints)
A2SegInd = findInterval(147782,IntervalPoints)
MPGSegInd = findInterval(65723,IntervalPoints)
POLSegInd = findInterval(29779,IntervalPoints)

#--------------------------
CF_1.Total.Neighbor.Count = GetNeighborCount(CF_1.Matrix)
# CF_1.Ratio = CF_1.Total.Neighbor.Count/sum(ContactFrequencyData_1)
CF_1.Ratio = sum(CF_1.Matrix[HS40SegInd,c(MPGSegInd,POLSegInd)])
cat ("--------------------------\n")
cat ("GM12878\n")
cat ("NeighborCount = ",CF_1.Total.Neighbor.Count,"\n")
cat ("max =", max(CF_1.Matrix),"\n")
cat ("sum =", sum(ContactFrequencyData_1),"\n")
cat ("ratio =", CF_1.Ratio,"\n")

#--------------------------
CF_2.Total.Neighbor.Count = GetNeighborCount(CF_2.Matrix)
# CF_2.Ratio = CF_2.Total.Neighbor.Count/sum(ContactFrequencyData_2)
CF_2.Ratio = sum(CF_2.Matrix[HS40SegInd,c(MPGSegInd,POLSegInd)])
cat ("K562\n")
cat ("NeighborCount = ",CF_2.Total.Neighbor.Count,"\n")
cat ("max =", max(CF_2.Matrix),"\n")
cat ("sum =", sum(ContactFrequencyData_2),"\n")
cat ("ratio =", CF_2.Ratio,"\n")

#----------------------------------------------------------------------------
CF_2.Matrix = CF_2.Matrix * (CF_1.Ratio/CF_2.Ratio)
max(CF_1.Matrix[A2SegInd,])
max(CF_2.Matrix[A2SegInd,])
max(CF_1.Matrix[HS40SegInd,])
max(CF_2.Matrix[HS40SegInd,])
#----------------------------------------
# histogram of A2 in GM and K cell
#----------------------------------------
LinkStr = ""
for (i in 1:(length(IntervalPoints)-1)){
	iRow = i
	PositionStr1.1 = IntervalPoints[iRow]
	PositionStr1.2 = IntervalPoints[iRow+1]
	Field1 = "hs16"
	params = sprintf("%.2f",CF_1.Matrix[i,A2SegInd])
	LinkStr = paste(LinkStr,paste(Field1," ",PositionStr1.1," ",PositionStr1.2," ",params,sep=""),sep="\n")
}

of_circos = paste(DestDir,"contactfreq.GM12878.A2Seg.Histgram.txt",sep="")
write.table(file=of_circos,LinkStr,quote=FALSE,row.names=FALSE,col.names=FALSE)

LinkStr = ""
for (i in 1:(length(IntervalPoints)-1)){
	iRow = i
	PositionStr1.1 = IntervalPoints[iRow]
	PositionStr1.2 = IntervalPoints[iRow+1]
	Field1 = "hs16"
	params = sprintf("%.2f",CF_2.Matrix[i,A2SegInd])
	LinkStr = paste(LinkStr,paste(Field1," ",PositionStr1.1," ",PositionStr1.2," ",params,sep=""),sep="\n")
}

of_circos = paste(DestDir,"contactfreq.K562.A2Seg.Histgram.txt",sep="")
write.table(file=of_circos,LinkStr,quote=FALSE,row.names=FALSE,col.names=FALSE)
#----------------------------------------
# histogram of HS40 in GM and K cell
#----------------------------------------
LinkStr = ""
for (i in 1:(length(IntervalPoints)-1)){
	iRow = i
	PositionStr1.1 = IntervalPoints[iRow]
	PositionStr1.2 = IntervalPoints[iRow+1]
	Field1 = "hs16"
	params = sprintf("%.2f",CF_1.Matrix[i,HS40SegInd])
	LinkStr = paste(LinkStr,paste(Field1," ",PositionStr1.1," ",PositionStr1.2," ",params,sep=""),sep="\n")
}

of_circos = paste(DestDir,"contactfreq.GM12878.HS40Seg.Histgram.txt",sep="")
write.table(file=of_circos,LinkStr,quote=FALSE,row.names=FALSE,col.names=FALSE)

LinkStr = ""
for (i in 1:(length(IntervalPoints)-1)){
	iRow = i
	PositionStr1.1 = IntervalPoints[iRow]
	PositionStr1.2 = IntervalPoints[iRow+1]
	Field1 = "hs16"
	params = sprintf("%.2f",CF_2.Matrix[i,HS40SegInd])
	LinkStr = paste(LinkStr,paste(Field1," ",PositionStr1.1," ",PositionStr1.2," ",params,sep=""),sep="\n")
}

of_circos = paste(DestDir,"contactfreq.K562.HS40Seg.Histgram.txt",sep="")
write.table(file=of_circos,LinkStr,quote=FALSE,row.names=FALSE,col.names=FALSE)


#----------------------------------------------------------------------------
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
if (CELL == "K562"){
	CF.Matrix = CF.Matrix * (CF_1.Ratio/CF_2.Ratio)
}
MaxCount = max(CF_1.Matrix,CF_2.Matrix)
CF.Matrix.Scale = CF.Matrix/MaxCount

#------------------------------------------------------------------
# Find interval starting index and end index of same sapce
#------------------------------------------------------------------
# GetStartEnd = function(SegPoints){
# 	NumPoints = length(SegPoints)
# 	NumSeg = NumPoints - 1
# 	MinSpace = min(SegPoints[2:NumPoints] - SegPoints[1:(NumPoints-1)])
# 	EquSpacePoints = seq(1,max(SegPoints),by=MinSpace)
#   PointsInd = findInterval(EquSpacePoints,SegPoints)
# 	StartEndInd.Vec = matrix(0,nrow=NumSeg,ncol=2)
# 	for (i in 1:NumSeg){
# 		StartEndInd.Vec[i,] = c(min(which(PointsInd==i)),max(which(PointsInd==i)))
# 	}
# 
# 	ConStartEndInd = cbind(1:nrow(StartEndInd.Vec),StartEndInd.Vec)
# 	ConStartEndSeg = cbind(EquSpacePoints, EquSpacePoints+MinSpace-1)
# 	return(list(ConStartEndInd=ConStartEndInd,ConStartEndSeg=ConStartEndSeg))
# }
# ConStartEnd = GetStartEnd(IntervalPoints)

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
RF.MaxCount = 10000

RF.Vec.Scale = RF.Vec = RFData
RF.Vec.Scale[,3] = RF.Vec.Scale[,3]/ RF.MaxCount

#------------------------------------------------------------------
# Ensemble contact for circos
#------------------------------------------------------------------

iCount = 1
LinkStr = ""
for (i in 1:nrow(RF.Vec)){
	iRow = RF.Vec[i,1]
	iCol = RF.Vec[i,2]
	if(RF.Vec[i,3]<0.001){next}
	PositionStr1.1 = IntervalPoints[iRow]
	PositionStr1.2 = IntervalPoints[iRow+1]
	PositionStr2.1 = IntervalPoints[iCol]
	PositionStr2.2 = IntervalPoints[iCol+1]
	Field1 = paste("segdup",sprintf("%04d hs16", iCount),sep="")
	freq = sprintf("thickness=%.3f",RF.Vec[i,3])
	maxfreq = sprintf("id=%.2f",RF.MaxCount,sep="")
	params = paste(freq,",",maxfreq,sep="")
	LinkStr = paste(LinkStr, paste(Field1," ",PositionStr2.1," ",PositionStr2.2," ",params,sep=""),sep="\n")
	LinkStr = paste(LinkStr, paste(Field1," ",PositionStr1.1," ",PositionStr1.2," ",params,sep=""),sep="\n")
	iCount = iCount + 1
}

of_circos = paste(DestDir,"contactfreq.ensemble.",SampleSize,".txt",sep="")
write.table(file=of_circos,LinkStr,quote=FALSE,row.names=FALSE,col.names=FALSE)

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

#----------------------------------------
# Propensity
#----------------------------------------

iCount = 1
LinkStr = ""
for (i in 1:nrow(ConInd_Vec)){
	iRow = ConInd_Vec[i,1]
	iCol = ConInd_Vec[i,2]
	if(PS.Vec[i,3]<0.001){next}
	PositionStr1.1 = IntervalPoints[iRow]
	PositionStr1.2 = IntervalPoints[iRow+1]
	PositionStr2.1 = IntervalPoints[iCol]
	PositionStr2.2 = IntervalPoints[iCol+1]
	Field1 = paste("segdup",sprintf("%04d hs16", iCount),sep="")
	freq = sprintf("thickness=%.3f",PS.Vec[i,3])
	maxfreq = sprintf("id=%.2f",PS.Max,sep="")
	params = paste(freq,",",maxfreq,sep="")
	LinkStr = paste(LinkStr, paste(Field1," ",PositionStr2.1," ",PositionStr2.2," ",params,sep=""),sep="\n")
	LinkStr = paste(LinkStr, paste(Field1," ",PositionStr1.1," ",PositionStr1.2," ",params,sep=""),sep="\n")
	iCount = iCount + 1
}

of_circos = paste(DestDir,"prop.",CELL,".",SampleSize,".txt",sep="")
write.table(file=of_circos,LinkStr,quote=FALSE,row.names=FALSE,col.names=FALSE)

#----------------------------------------
# Linker of A2 in GM and K cell
#----------------------------------------
iCount = 1
LinkStr = ""
A2ConInd = which(PS.Vec[,2]==A2SegInd | PS.Vec[,1]==A2SegInd)
for (i in 1:(length(A2ConInd))){
	iRow = PS.Vec[A2ConInd[i],1]
	iCol = PS.Vec[A2ConInd[i],2]
	PositionStr1.1 = IntervalPoints[iRow]
	PositionStr1.2 = IntervalPoints[iRow+1]
	PositionStr2.1 = IntervalPoints[iCol]
	PositionStr2.2 = IntervalPoints[iCol+1]
	Field1 = paste("segdup",sprintf("%04d hs16", iCount),sep="")
	freq = sprintf("thickness=%.3f",PS.Vec[A2ConInd[i],3])
	maxfreq = sprintf("id=%.2f",PS.Max,sep="")
	params = paste(freq,",",maxfreq,sep="")
	LinkStr = paste(LinkStr, paste(Field1," ",PositionStr2.1," ",PositionStr2.2," ",params,sep=""),sep="\n")
	LinkStr = paste(LinkStr, paste(Field1," ",PositionStr1.1," ",PositionStr1.2," ",params,sep=""),sep="\n")
	iCount = iCount + 1
}
of_circos = paste(DestDir,"contactfreq.",CELL,".A2Seg.Con.Link.txt",sep="")
write.table(file=of_circos,LinkStr,quote=FALSE,row.names=FALSE,col.names=FALSE)

#----------------------------------------
# Linker of HS40 in GM and K cell
#----------------------------------------
iCount = 1
LinkStr = ""
HS40ConInd = which(PS.Vec[,2]==HS40SegInd | PS.Vec[,1]==HS40SegInd)
for (i in 1:(length(HS40ConInd))){
	iRow = PS.Vec[HS40ConInd[i],1]
	iCol = PS.Vec[HS40ConInd[i],2]
	PositionStr1.1 = IntervalPoints[iRow]
	PositionStr1.2 = IntervalPoints[iRow+1]
	PositionStr2.1 = IntervalPoints[iCol]
	PositionStr2.2 = IntervalPoints[iCol+1]
	Field1 = paste("segdup",sprintf("%04d hs16", iCount),sep="")
	freq = sprintf("thickness=%.3f",PS.Vec[HS40ConInd[i],3])
	maxfreq = sprintf("id=%.2f",PS.Max,sep="")
	params = paste(freq,",",maxfreq,sep="")
	LinkStr = paste(LinkStr, paste(Field1," ",PositionStr2.1," ",PositionStr2.2," ",params,sep=""),sep="\n")
	LinkStr = paste(LinkStr, paste(Field1," ",PositionStr1.1," ",PositionStr1.2," ",params,sep=""),sep="\n")
	iCount = iCount + 1
}
of_circos = paste(DestDir,"contactfreq.",CELL,".HS40Seg.Con.Link.txt",sep="")
write.table(file=of_circos,LinkStr,quote=FALSE,row.names=FALSE,col.names=FALSE)






#------------------------------------------------------------------
# Bootstrap for p-value
#------------------------------------------------------------------
# rebuild the distribution
BootstrapTest = function(xcount,xtotal,ycount,ytotal,alpha){
	# i=2
	# xcount = CF.Matrix[ConInd_Vec[i,1],ConInd_Vec[i,2]]
	# xtotal = MaxCount
	# ycount = RF.Vec[i,3]
	# ytotal = RF.MaxCount	
	
	if (ycount ==0){ycount = 1}
	x_mean = xcount/xtotal
	y_mean = ycount/ytotal
	sampleid = sample(ytotal,ycount)
	ReadCount = rep(0,ytotal)
	ReadCount[sampleid] = 1
	
	N_Trial = 30
	Mean_Vec = rep(0,N_Trial)
	for(k in 1:N_Trial){
		# cat("\t",k,"/",N_Trial,"\n")
		BootStrapIDs =sample(ytotal,replace=T)
		Mean_Vec[k] = mean(ReadCount[BootStrapIDs])
	}
	Mean_Boot = mean(Mean_Vec)
	SE_Boot = sqrt((1/(N_Trial-1))* sum((Mean_Vec-Mean_Boot)^2))
	t_star = qt(c(alpha,1-alpha),df=N_Trial-1)
	CI_Boot = Mean_Boot+t_star*SE_Boot
	# cat(CI_Boot,"\n")
	
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
# cat("Calculating real P-value...\n")
PValue.Vec = cbind(ConInd_Vec,rep(0,nrow(ConInd_Vec)),rep(0,nrow(ConInd_Vec)))
for (i in 1:nrow(ConInd_Vec)){
	# cat(i,"/",nrow(ConInd_Vec),"\n")
	# cat(RF.Vec[i,3],"\n")
 	BootStrapTest1 = BootstrapTest(CF.Matrix[ConInd_Vec[i,1],ConInd_Vec[i,2]],MaxCount, RF.Vec[i,3], RF.MaxCount,alpha)
	PValue.Vec[i,3] = BootStrapTest1$P_Value
	PValue.Vec[i,4] = BootStrapTest1$l_e_g
}

# pvalueDir = "data/analysis/ENm008/bootstrap/"
# of_txt = paste(pvalueDir,"bootstrap_",CELL,"_pvalue.txt",sep="")
# # write.table(file=of_txt,PValue.Vec,quote=FALSE,row.names=FALSE,col.names=FALSE)
# PValue.Vec = read.table(of_txt)


#--------------------
# FDR (False Discovery Rate) by cluster calculation
#--------------------
# cat("Calculating FDR...\n")

N_FDR = 30
fdr_Matrix = matrix(0,nrow=nrow(ConInd_Vec),ncol=N_FDR)
for (i_fdr in 1:N_FDR){
	ShuffeledIDs = sample(nrow(ConInd_Vec))
	# cat("fileid =",fileid,"i_fdr = ", i_fdr,"\n")
	for (i_con in 1:nrow(ConInd_Vec)){
   # cat("\t i_con =",i_con,"/", "totalcon =",nrow(ConInd_Vec),"\n")
		BootStrapTest2 = BootstrapTest(CF.Matrix[ConInd_Vec[ShuffeledIDs[i_con],1],ConInd_Vec[ShuffeledIDs[i_con],2]],MaxCount, RF.Vec[i_con,3], RF.MaxCount,alpha)
		fdr_Matrix[i_con,i_fdr] = BootStrapTest2$P_Value
	}
}

# fdr_path = "data/analysis/ENm008/fdr/"
# of_txt = paste(fdr_path,"fdr_",CELL,"_",fileid,".txt",sep="")
# write.table(file=of_txt,formatC(fdr_Matrix,format="e"),quote=F,col.names=F,row.names=F)

#--------------------
# Read FDR files, and get significant
#--------------------
# 
# fdr_path = "data/analysis/ENm008/fdr/"
# fdr_files = list.files(path=fdr_path,pattern=paste("fdr_",CELL,sep=""))
# fdr_Matrix = read.table(paste(fdr_path,fdr_files[1],sep=""))
# for (i in 2:length(fdr_files)){
# 	fdr_Matrix = cbind(fdr_Matrix,read.table(paste(fdr_path,fdr_files[i],sep="")))
# }
# 

V_p = apply(fdr_Matrix,2,function(x){length(which(x<alpha))})
s_p = length(which(PValue.Vec[,3]<alpha))
fdr = mean(V_p/(V_p+s_p))
PValue.Vec.Sort.Ind = sort.int(PValue.Vec[,3],index.return=T)$ix
PickSize = round(s_p * (1-fdr))

cat("CELL=",CELL,"\n")
cat("Total Contact=",nrow(ConInd_Vec),"\n")
cat("Original Sig Cont =", s_p,"\n")
cat("Median of PS",median(PS.Vec[,3]),"\n")
cat("FDR =",fdr,"\n")
cat("PickSize =",PickSize,"\n")


#--------------------
# Output to circos
#--------------------

iCount = 1
LinkStr = ""
for (k in 1:PickSize){
	i = PValue.Vec.Sort.Ind[k]
	iRow = PValue.Vec[i,1]
	iCol = PValue.Vec[i,2]
	if(PValue.Vec[i,4] != 2){next}
	PositionStr1.1 = IntervalPoints[iRow]
	PositionStr1.2 = IntervalPoints[iRow+1]
	PositionStr2.1 = IntervalPoints[iCol]
	PositionStr2.2 = IntervalPoints[iCol+1]
	Field1 = paste("segdup",sprintf("%04d hs16", iCount),sep="")
	freq = sprintf("thickness=%.3f",PS.Vec[i,3])
	maxfreq = sprintf("id=%.2f",PS.Max,sep="")
	params = paste(freq,",",maxfreq,sep="")
	LinkStr = paste(LinkStr, paste(Field1," ",PositionStr2.1," ",PositionStr2.2," ",params,sep=""),sep="\n")
	LinkStr = paste(LinkStr, paste(Field1," ",PositionStr1.1," ",PositionStr1.2," ",params,sep=""),sep="\n")
	iCount = iCount + 1
}

of_circos = paste(DestDir,"pval.pos.",CELL,".",alpha,".",SampleSize,".txt",sep="")
write.table(file=of_circos,LinkStr,quote=FALSE,row.names=FALSE,col.names=FALSE)
# 
iCount = 1
LinkStr = ""
for (k in 1:PickSize){
	i = PValue.Vec.Sort.Ind[k]
	iRow = PValue.Vec[i,1]
	iCol = PValue.Vec[i,2]
	if(PValue.Vec[i,4] != 0){next}
	PositionStr1.1 = IntervalPoints[iRow]
	PositionStr1.2 = IntervalPoints[iRow+1]
	PositionStr2.1 = IntervalPoints[iCol]
	PositionStr2.2 = IntervalPoints[iCol+1]
	Field1 = paste("segdup",sprintf("%04d hs16", iCount),sep="")
	freq = sprintf("thickness=%.3f",PS.Vec[i,3])
	maxfreq = sprintf("id=%.2f",PS.Max,sep="")
	params = paste(freq,",",maxfreq,sep="")
	LinkStr = paste(LinkStr, paste(Field1," ",PositionStr2.1," ",PositionStr2.2," ",params,sep=""),sep="\n")
	LinkStr = paste(LinkStr, paste(Field1," ",PositionStr1.1," ",PositionStr1.2," ",params,sep=""),sep="\n")
	iCount = iCount + 1
}

of_circos = paste(DestDir,"pval.neg.",CELL,".",alpha,".",SampleSize,".txt",sep="")
write.table(file=of_circos,LinkStr,quote=FALSE,row.names=FALSE,col.names=FALSE)

#--------------------
# Output to calculation
#--------------------
PVal.Output = PValue.Vec[PValue.Vec.Sort.Ind[1:PickSize],]
# sort first segment
PVal.Output = PVal.Output[sort.int(PVal.Output[,1],index.return=T)$ix,]
# sort second segment
PVal.Output = PVal.Output[sort.int(PVal.Output[,2],index.return=T)$ix,]
# shorten output of pval
PVal.Output[,3] = format(PVal.Output[,3],digits=3,scientific=T)
# file name
of_name = paste(DestDir,CELL,".",SampleSize,".pval.",alpha,".txt",sep="")
# output segind, segind, pvalue and attraction(2) or repulsion(0)
# write.table(file=of_name,
# 	t(c(nrow(CF.Matrix),ncol(CF.Matrix))),quote=FALSE,row.names=FALSE,col.names=FALSE )
# write.table(file=of_name,
# 	 PVal.Output[,c(1,2,4)],quote=FALSE,row.names=FALSE,col.names=FALSE,append=T)
write.table(file=of_name,
	 PVal.Output[,c(1,2,4)],quote=FALSE,row.names=FALSE,col.names=FALSE)


