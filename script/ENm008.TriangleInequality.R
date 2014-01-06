#!/usr/bin/R
rm(list=ls())
SegmentFN = "data/analysis/ENm008/ENm008.p1500.equ.seg.len.txt"
# ContactFN = "data/analysis/ENm008/GM12878.p1500_d3300.q5.a5.sis.node.dst"
ContactFN = "data/analysis/ENm008/K562.p1500_d3300.q5.a5.sis.node.dst"

SegLenVec = read.table(SegmentFN)
ContactVec = read.table(ContactFN)

ContactVec = rbind(ContactVec , cbind(1:nrow(SegLenVec),2:(nrow(SegLenVec)+1),SegLenVec[,1]))
# sort

ContactVec = ContactVec[sort.int(ContactVec[,1],index.return=T)$ix,]
ContactVec = ContactVec[sort.int(ContactVec[,2],index.return=T)$ix,]
ContactVec = as.matrix(ContactVec)
row.names(ContactVec) = NULL

findTriangleInd <- function(Segment){
	# Segment = ContactVec[82,]
	A_Ind = Segment[1]
	B_Ind = Segment[2]
	AB_Len = Segment[3]
	AorBInd = which(ContactVec[,1]==A_Ind | ContactVec[,1]== B_Ind)
	RowIndVec = AorBInd[which(duplicated(ContactVec[AorBInd,2]))]
	C_IndVec = ContactVec[RowIndVec,2]
	
	ConInd = NULL
	for (i in 1:length(C_IndVec)){
		C_Ind = C_IndVec[i]
		ConInd = rbind(ConInd,which((ContactVec[,1]==A_Ind | ContactVec[,1]== B_Ind) & ContactVec[,2]==C_Ind))
	}
	return(ConInd)
}

TriangleIndList = apply(ContactVec,1,
	function(x){
		CIndVec=findTriangleInd(x); 
		if(length(CIndVec)>0){
			return(CIndVec)
		}
		else{
			return(0)
		}})
		
TriangleIndList [[88]]

for ( i in 1:length(TriangleIndList)){
	iNumRow = nrow(TriangleIndList[[i]])
	if (is.null(iNumRow)){ next }
	for (iRow in 1:iNumRow){
		Ind_1 = i
		Ind_2 = TriangleIndList[[i]][iRow,1]
		Ind_3 = TriangleIndList[[i]][iRow,2]
		AB_Len = ContactVec[Ind_1,3] 
		AC_Len = ContactVec[Ind_2,3] 
		BC_Len = ContactVec[Ind_3,3] 
		# if (abs(AB_Len - AC_Len)> BC_Len - 300 | abs(AB_Len - BC_Len)>AC_Len -300| abs(BC_Len - AC_Len)>AB_Len-300){
			if (abs(AB_Len - AC_Len)> BC_Len | abs(AB_Len - BC_Len)>AC_Len| abs(BC_Len - AC_Len)>AB_Len){
			cat ("------\n",i,"\n")
			# cat (ContactVec[Ind_1,] ,"\n",ContactVec[Ind_2,] ,"\n",ContactVec[Ind_3,] ,"\n",sep=' ')
		}
	}
}

