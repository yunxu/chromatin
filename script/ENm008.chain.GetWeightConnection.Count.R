#!/usr/bin/R
rm(list=ls())
#------------------------------------------------------------------
library(grid)
library(gridExtra)
library(lattice)
library(RColorBrewer)
library(reshape)
library(ggplot2)
#------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (Sys.info()["sysname"] == "Darwin"){
	# CELL = "GM12878"
	CELL = "K562"
}else{
	CELL = args[1]
}

#------------------------------------------------------------------
# NumMax = 1000
Type="r270.ovl"
method="alpha"
#------------------------------------------------------------------
Type="con"
method = "cut.all"
#------------------------------------------------------------------
ConDir=paste("/Users/yunxu/workspace/projects/chromatin/result/analysis/ENm008/chain/",method,"/par.1_1_1_1/",CELL,".p1500_d3300",sep="")
# ConDir=paste("/dump/yunxu/working/projects.new/chromatin/result/analysis/ENm008/chain/",method,"/par.1_1_1_1/",CELL,".p1500_d3300",sep="")
#------------------------------------------------------------------
GetWeightFile <- function(ConDir){
	ConFiles   = list.files(path=ConDir,pattern=Type,full.names=T)
	NumFiles = length(ConFiles)
# NumFiles = 10
	Weight = rep(0,NumFiles)
	ConList = list()
	for (i in 1:NumFiles){
		Weight[i] = read.table(ConFiles[i],comment.char="",nrows=1)[,3]
		ConList[[i]] = read.table(ConFiles[i])
	}
	return(list("weight"=Weight,"conlist"=ConList))
}
#------------------------------------------------------------------
SegLenFN="data/analysis/ENm008/ENm008.p1500.equ.seg.len.txt"
NumNodes=nrow(read.table(SegLenFN))+1
#------------------------------------------------------------------

# WeightData <- GetWeightFile(ConDir)
# Weight <- WeightData$weight
# ConList <- WeightData$conlist
# NumFiles <- length(Weight)
#------------------------------------------------------------------
# NewPIJMat = matrix(NA,nrow=NumNodes,ncol=NumNodes)
# for (i in 1:(NumNodes-2)){
# 	for (j in (i+2):NumNodes){
# 		cat(i, " / ",j ," / ", NumNodes,"\n"); flush.console();
# 		ConIJ = unlist(lapply(ConList,function(x){length(which(x[,1]==i & x[,2]==j))}))
# 		NewPIJMat[i,j] = sum(ConIJ * Weight) 
# 	}
# }
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
	}
	# M = NewPIJMat
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

# MIN_M = 0
# MAX_M = 1000
# # NewPIJMat_Scale = NewPIJMat *1000/max(NewPIJMat,na.rm=T)
# # myPalette <- brewer.pal(11, "Spectral")
# myPalette <- brewer.pal(9,"YlGnBu")
# # FN = paste(OutDir,"/","w_",NumNodesofOneWindow,"_mean.pdf",sep="")
# FN = paste(OutDir,"/",CELL,".p1500_d3300.c840.pij.raw.count.pdf",sep="")
# cat(FN,"\n")
# DrawTriangleHeatmap(FN,NewPIJMat_Scale,MIN_M,MAX_M,myPalette,"All")


read.triplet.matrix <- function(file){
	# file = FN
	HeadInfo = read.table(file,comment.char="",nrows=1)
	NumRow = HeadInfo[1,2]
	NumCol = HeadInfo[1,3]

	m = matrix(0,nrow=NumRow,ncol=NumCol)
	DF = read.table(file)
	for (i in 1:nrow(DF)){
		ind1 = which(DF[i,1] == Equ_Contact_Point_Ind[,1])
		ind2 = which(DF[i,2] == Equ_Contact_Point_Ind[,2])
		matched = as.logical(length(intersect(ind1,ind2)))
		if (matched){
			# m[DF[i,1],DF[i,2]] = DF[i,3]
			# m[DF[i,2],DF[i,1]] = DF[i,3]
			# if(DF[i,1] <20 & DF[i,2] >30 & DF[i,2]<50){
			# 	m[DF[i,1],DF[i,2]] = m[DF[i,2],DF[i,1]] = 0
			# }
			m[DF[i,1],DF[i,2]] = NA
			m[DF[i,2],DF[i,1]] = NA
		} else {
			m[DF[i,1],DF[i,2]] = DF[i,3]
			m[DF[i,2],DF[i,1]] = DF[i,3]
			# m[DF[i,1],DF[i,2]] = NA
			# m[DF[i,2],DF[i,1]] = NA
		}
	}

	return(m)
}

(FN <- paste("tmp/",CELL,".p1500_d3300.c840.",Type,".pij.triplet.allall.txt",sep=""))
mat <- read.triplet.matrix(FN)
mat <- mat * 1000/max(mat,na.rm=T)
FN = paste(OutDir,"/",CELL,".p1500_d3300.c840.con.pij.allall.pdf",sep="")
FN = paste(OutDir,"/",CELL,".p1500_d3300.c840.con.pij.common.pdf",sep="")
FN = paste(OutDir,"/",CELL,".p1500_d3300.c840.con.pij.new.pdf",sep="")
cat(FN,"\n")
MIN_M = 0
MAX_M = 1000
myPalette <- brewer.pal(9,"YlGnBu")
myPalette <- brewer.pal(9,"YlOrBr")
myPalette <- brewer.pal(9,"PuBu")
myPalette[1:2] = myPalette[3]
DrawTriangleHeatmap(FN,mat,MIN_M,MAX_M,myPalette,"All")
