#!/usr/bin/R

rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
# Par = "par.1_1_1_1"
# SampleSize = "GM12878.p1500_d3300"
# SampleSize = "K562.p1500_d3300"
# 
# NumNodesofOneWindow = 1
# MaxNum=1000

 Par = args[1]
 SampleSize = args[2]
 NumNodesofOneWindow = as.numeric(args[3])
 MaxNum = as.numeric(args[4])

CELL = gsub("(\\w*).(\\w*)","\\1",SampleSize)
DataDir = paste(Par,"/",SampleSize,sep="")
OutDir = paste("/tmp/",DataDir,"/",sep="")

library(grid)
library(lattice)
library(RColorBrewer)

FileNames = list.files(DataDir,pattern="*.pts",recursive=T,full.names=T)
#FileNames = FileNames[1:5]

LogWeight = rep(0,length(FileNames))
for (i in 1:length(FileNames)){
	FileName = FileNames[i]
  LogWeight[i] = read.table(FileName,comment.char="",nrows=1)[,3]
}

MaxLogWeight = max(LogWeight)
Weight = exp(LogWeight - MaxLogWeight)
NumCopies = round (Weight * MaxNum)

NonZeroInd = which(NumCopies!=0)
NumCopyVec = NumCopies[NonZeroInd]
NewFileNames = FileNames[NonZeroInd]

PtsSample = read.table(NewFileNames[1])
NumNodes = nrow(PtsSample)
PtsArr = array(0, dim=c(NumNodes,3, length(NewFileNames)))

for ( i in 1:length(NewFileNames)){
	# cat(i,"\n")
	flush.console()
	FN = NewFileNames[i]
	PtsArr[,,i] = as.matrix(read.table(FN)[,1:3])
}
PtsArr = PtsArr / 10

NumWindows = NumNodes - NumNodesofOneWindow + 1
WinPtsArr = array(0, dim=c(NumWindows,3,length(NewFileNames)))
for (i in 1:length(NewFileNames)){
	# TmpPtsMat = matrix(0,nrow=NumWindows,ncol=3)
	for ( j in 1:NumWindows){
		StartInd = j
		EndInd   = j + NumNodesofOneWindow - 1
		WinPtsArr[j,,i] = apply(PtsArr[StartInd:EndInd,,i,drop=F],2,mean)
	}
}



NumNodeNodeCon = sum(1:(NumWindows-1))
# DistMat = matrix(0,nrow=NumNodeNodeCon,ncol=length(NewFileNames))
NodeNodeVec = matrix(0,nrow=NumNodeNodeCon,ncol=2)
DistArr     = array(0, dim=c(NumWindows,NumWindows,length(NewFileNames)))
MedianMat   = matrix(NA,nrow=NumWindows,ncol=NumWindows)
MeanMat   = matrix(NA,nrow=NumWindows,ncol=NumWindows)
# VarMat      = matrix(NA,nrow=NumWindows,ncol=NumWindows)
SDMat      = matrix(NA,nrow=NumWindows,ncol=NumWindows)
k           = 1
for (i in 2:NumWindows){
	cat(i,"/",NumWindows,"\n")
	for (j in  1:(i-1)){
		flush.console()
		NodeNodeVec[k,] = c(j,i)
		DistArr[j,i,] = apply(WinPtsArr[c(j,i),,],3,dist)
		
		MedianMat[j,i] = median(DistArr[j,i,])
		MeanMat[j,i] = mean(DistArr[j,i,])
		# VarMat[j,i] = var(DistArr[j,i,])
		SDMat[j,i] = sd(DistArr[j,i,])
		k = k + 1
	}
}


convertToColors <- function(mat, min.dist, max.dist, Palette) {
    # Produce 'normalized' version of matrix, with values ranging from 0 to 1
    rng <- range(mat, na.rm = TRUE)
		minM = min(min.dist,rng[1])
		maxM = max(max.dist,rng[2])
    m <- (mat - minM)/(maxM - minM)
    # Convert to a matrix of sRGB color strings
    m2 <- m; 
		class(m2) <- "character"
    # m2[!is.na(m2)] <- rgb(colorRamp(heat.colors(16))(m[!is.na(m)]), max = 255)
    m2[!is.na(m2)] <- rgb(colorRamp(Palette)(m[!is.na(m)]), max = 255)
    m2[is.na(m2)] <- "transparent"
    return(m2)
}
# http://colorbrewer2.org/
# display.brewer.all()

# myPalette <- rev(brewer.pal(9, "PuBuGn"))
# myPalette <- rev(brewer.pal(9, "BuGn"))
## Initialize plot and prepare two viewports
# grid.newpage()
Scale2Dist <- function(c) {round(MAX_M *c/max(c))}

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
								text.cex = 0.5
							)
	  },
	  axis.default(side = side, ...))
}

scaleMatrix <- function(mat, min_m, max_m){
	m <- (mat - min_m) / (max_m - min_m)
	return(m)
}


DrawTriangleHeatmap <- function(FN,Mat,MIN_M, MAX_M, Palette,type){

	pdf(FN)	
	# M = MedianMat
	# MIN_M = MIN_Dist
	# MAX_M = MAX_Dist
	# Palette = myPalette
	M = scaleMatrix(Mat,MIN_M, MAX_M)
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
				for (i in 1:nrow(M)){
					panel.text(i+0.25,i-0.5,i,cex=.5,srt=90)
				}
				panel.text(nrow(M)/2+3,nrow(M)/2-3,paste("Window size =",NumNodesofOneWindow),cex=1,srt=45)
				panel.text(nrow(M)/2+5,nrow(M)/2-5,type,cex=1,srt=45)
				panel.text(nrow(M)/2+7,nrow(M)/2-7,CELL,cex=1,srt=45)
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
	legend.vp <- viewport(	x = unit(2.5,"lines"),	y = unit(0.7,"npc"),
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

MIN_M = 0
MAX_M = 400
myPalette <- brewer.pal(11, "Spectral")
FN = paste(OutDir,"/","w_",NumNodesofOneWindow,"_mean.pdf",sep="")
cat(FN,"\n")
DrawTriangleHeatmap(FN,MeanMat,MIN_M,MAX_M,myPalette,"Mean")

MIN_M = 0
MAX_M = 400
myPalette <- brewer.pal(11, "Spectral")
FN = paste(OutDir,"/","w_",NumNodesofOneWindow,"_median.pdf",sep="")
cat(FN,"\n")
DrawTriangleHeatmap(FN,MedianMat,MIN_M,MAX_M,myPalette,"Median")

MIN_M = 0
MAX_M = 160
myPalette <- brewer.pal(11, "RdYlBu")
FN = paste(OutDir,"/","w_",NumNodesofOneWindow,"_sd.pdf",sep="")
cat(FN,"\n")
DrawTriangleHeatmap(FN,SDMat,MIN_M,MAX_M,myPalette,"Standard deviation")
