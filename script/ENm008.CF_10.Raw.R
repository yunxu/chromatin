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

Equ_Contact_Point_Ind = cbind(
  GetEquPointIndFromExpPointPos(CFData_1[,1]),
  GetEquPointIndFromExpPointPos(CFData_1[,2]))
  
# remove equ neighbor segment
EquSegNonNeighborInd = which(Equ_Contact_Point_Ind[,2]-Equ_Contact_Point_Ind[,1]!=1)
Equ_Contact_Point_Ind = Equ_Contact_Point_Ind[EquSegNonNeighborInd,]

CFData_1 = CFData_1[EquSegNonNeighborInd,]
CFData_2 = CFData_2[EquSegNonNeighborInd,]

FN_1 = "/tmp/GM12878.Raw.Delete.Short.Connection.Circos.txt"
WriteCircos_RawConnection(FN_1, CFData_1)
FN_2 = "/tmp/K562.Raw.Delete.Short.Connection.Circos.txt"
WriteCircos_RawConnection(FN_2, CFData_2)


CF_1.Vec = cbind(Equ_Contact_Point_Ind,CFData_1[,3])
CF_2.Vec = cbind(Equ_Contact_Point_Ind,CFData_2[,3])

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
#----------------------------------------------------------------------------

# CFData_1[CFData_1[,1]==167103 | CFData_1[,2]==167103, ]
# CFData_2[CFData_2[,1]==167103 | CFData_2[,2]==167103, ]

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


#------------------------------------------------------------------
# Propensity
#------------------------------------------------------------------
PS_1.Vec = cbind(Con.Vec,
	apply(RF.Vec.Scale[,3:ncol(RF.Vec.Scale)], 2, function(x){CF_1.Vec.Scale[,3]/x}))
PS_2.Vec = cbind(Con.Vec,
	apply(RF.Vec.Scale[,3:ncol(RF.Vec.Scale)], 2, function(x){CF_2.Vec.Scale[,3]/x}))











#------------------------------------------------------------------
# Draw triangle heatmap
#------------------------------------------------------------------
library(grid)
library(gridExtra)
library(lattice)
library(RColorBrewer)
library(reshape)
library(ggplot2)
#------------------------------------------------------------------
scaleMatrix <- function(mat, min_m, max_m){
	m <- (mat - min_m) / (max_m - min_m)
	return(m)
}
#------------------------------------------------------------------
Scale2Dist <- function(c) {sprintf("%.0f",MAX_M *c/max(c))}
#------------------------------------------------------------------

axis.legend <- function(side, ...) {
	ylim <- current.panel.limits()$ylim
	switch(side,
	  left = {
			prettyY <- pretty(ylim)
			labY <- prettyY
			panel.axis(side = side, outside = TRUE,
			           at = prettyY, labels = labY,
								draw.labels =F,ticks=F)
	  },
	  right = {
			prettyY <- pretty(ylim)
			# seqY <- seq(1,max(prettyY))
			seqY <- prettyY
			labY <- Scale2Dist(seqY)
			
			panel.axis(side = side, outside = T,
			           at = seqY,
									# line.col = "black",
								 labels = labY,
								text.cex = 1.5
							)
	  },
	  axis.default(side = side, ...))
}

#------------------------------------------------------------------

DrawTriangleHeatmap <- function(FN,Mat,MIN_M, MAX_M, Palette,type){

	library(tools)
	FileExt = file_ext(FN)
	if (FileExt == "ps"){
		postscript(FN)
	} else {
		pdf(FN)	
	}	# M = NewPIJMat
	# MIN_M = MIN_Dist
	# MAX_M = MAX_Dist
	# Palette = myPalette
	M = scaleMatrix(Mat,MIN_M, MAX_M)
	M[which(lower.tri(M))] = NA
	# trellis.par.set(axis.line=list(lwd=0,col="white"))
	HeatMapPlot <- levelplot(
		M,
		at = do.breaks(range(0,1),length(Palette)),
		col.regions=colorRampPalette(Palette, space = "Lab"),
		colorkey= F,
		xlab = "", 
		ylab="",
		scales = list(draw = FALSE),
		par.settings = list(axis.line=list(lwd=0,col="white")),
		panel = function(...){
		    panel.levelplot(...)
				for (i in (1:as.integer(nrow(M)/10))*10){
					panel.text(i+0.25,i-0.5,i,cex=2,srt=90)
				}
				# panel.text(nrow(M)/2+3,nrow(M)/2-3,paste("Window size =",NumNodesofOneWindow),cex=1,srt=45)
				# panel.text(nrow(M)/2+5,nrow(M)/2-5,type,cex=1,srt=45)
				# panel.text(nrow(M)/2+7,nrow(M)/2-7,CELL,cex=1,srt=45)

				step = 5
				lineX = seq(step,nrow(M),by=step)
				for (i in 1:length(lineX)){
					panel.lines(
						x=c(lineX[i],lineX[i]),
						y=c(lineX[i],nrow(M)),
						col=grey(0.9),lwd=0.5)
					panel.lines(
						x=c(1,lineX[i]),
						y=c(lineX[i],lineX[i]),
						col=grey(0.9),lwd=0.5)
				}
				
		}
	)
	
	# trellis.par.set(axis.line=list(lwd=1))
	LegendPlot <- levelplot(
		t(as.matrix(do.breaks(range(MIN_M,MAX_M),length(Palette)-1))),
		col.regions=colorRampPalette(Palette, space = "Lab"),
		colorkey=F,
		xlab = "", ylab="",
		scales = list(
			x = list(draw=F)),
			par.settings = list(axis.line=list(lwd=.25,col="black")),
		axis = axis.legend,
	)
	
	
	
	
	all.layout = grid.layout(3,3,
	  widths = unit(c(6,1,5), c("lines","null","lines")), 
	  heights = unit(c(6,1,5), c("lines", "null", "lines")))

	top.vp <- viewport(layout = all.layout)

	margin1.vp <- viewport(layout.pos.col = 2, layout.pos.row = 3, name = "margin1")
	margin2.vp <- viewport(layout.pos.col = 1, layout.pos.row = 2, name = "margin2")
	margin3.vp <- viewport(layout.pos.col = 2, layout.pos.row = 1, name = "margin3")
	margin4.vp <- viewport(layout.pos.col = 3, layout.pos.row = 2, name = "margin4")
	plot.vp <- viewport(layout.pos.col = 2, layout.pos.row = 2, name = "plot")

	heatmap.vp <- viewport(angle = -45, name="heatmap")
	legend.vp <- viewport(	x = unit(0.5,"lines"),	y = unit(0.9,"npc"),
		width = unit(.25,"lines"),
		height = unit(8,"lines"),
		name = "legend"
		)
	n <- nrow(M) 

	label.vp <- viewport(name="label",width=sqrt(2)*(n-1)/n ,
	 height=unit(1,"lines"))


	plot.vptree <- vpTree(plot.vp,vpList(heatmap.vp,label.vp))
	# top.vptree <- vpTree(top.vp, vpList(margin1.vp, margin2.vp, margin3.vp, margin4.vp, plot.vptree))
	top.vptree <- vpTree(top.vp,
		vpList(margin1.vp,margin2.vp,margin3.vp,margin4.vp,plot.vptree))

	grid.newpage()
	# grid.show.layout(all.layout)
	pushViewport(top.vptree)
	seekViewport("heatmap")
	# seekViewport("plot")
		 ow <- options("warn")
		 options(warn = -1)
		trellis.par.set(axis.line=list(lwd=0))
		print(HeatMapPlot,vp=heatmap.vp,newpage=FALSE)
		 options(ow)
	  # grid.raster(t(convertToColors(M,MIN_M,MAX_M,Palette)), interpolate = F)
	upViewport()

	seekViewport("margin4")
	  # grid.rect(gp = gpar(col = "gray90"))
		pushViewport( legend.vp )
		trellis.par.set(axis.line=list(lwd=1))
		print(LegendPlot,newpage=FALSE)
			# seekViewport("legend")		
	upViewport(0)
	dev.off()
}
#------------------------------------------------------------------
OutDir = "/tmp"
#------------------------------------------------------------------

#------------------------------------------------------------------
NumNodes = 54
CELL = "GM12878"
CF_Sample = CFData_1
CELL = "K562"
CF_Sample = CFData_2

NewPIJMat = matrix(NA,nrow=NumNodes,ncol=NumNodes)
# CF_Sample = CF_2.Vec
for (i in 1:nrow(CF_Sample)){
	# NewPIJMat[CF_Sample[i,1],CF_Sample[i,2]] = CF_Sample[i,3]
		NewPIJMat[Equ_Contact_Point_Ind[i,1],Equ_Contact_Point_Ind[i,2]] = CF_Sample[i,3]
}
#------------------------------------------------------------------
MIN_M = 0
MAX_M = 1000
# myPalette <- brewer.pal(11, "Spectral")
myPalette <- brewer.pal(9,"YlGnBu")
myPalette[1:2] = myPalette[3]
# FN = paste(OutDir,"/","w_",NumNodesofOneWindow,"_mean.pdf",sep="")
FN = paste(OutDir,"/",CELL,".Raw.Triangle.Heatmap.pdf",sep="")
cat(FN,"\n")
DrawTriangleHeatmap(FN,NewPIJMat,MIN_M,MAX_M,myPalette,"All")
