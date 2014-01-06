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
	CELL = "K562"
	# CELL = "GM12878"
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

WeightData <- GetWeightFile(ConDir)
Weight <- WeightData$weight
ConList <- WeightData$conlist
NumFiles <- length(Weight)
#------------------------------------------------------------------
NewPIJMat = matrix(NA,nrow=NumNodes,ncol=NumNodes)
for (i in 1:(NumNodes-2)){
	for (j in (i+2):NumNodes){
		cat(i, " / ",j ," / ", NumNodes,"\n"); flush.console();
		ConIJ = unlist(lapply(ConList,function(x){length(which(x[,1]==i & x[,2]==j))}))
		NewPIJMat[i,j] = sum(ConIJ * Weight) / sum(Weight)
	}
}
#------------------------------------------------------------------
scaleMatrix <- function(mat, min_m, max_m){
	m <- (mat - min_m) / (max_m - min_m)
	return(m)
}
#------------------------------------------------------------------
Scale2Dist <- function(c) {sprintf("%.1f",MAX_M *c/max(c))}
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
								text.cex = 0.5
							)
	  },
	  axis.default(side = side, ...))
}

#------------------------------------------------------------------

DrawTriangleHeatmap <- function(FN,Mat,MIN_M, MAX_M, Palette,type){

	pdf(FN)	
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
				for (i in 1:nrow(M)){
					panel.text(i+0.25,i-0.5,i,cex=.5,srt=90)
				}
				# panel.text(nrow(M)/2+3,nrow(M)/2-3,paste("Window size =",NumNodesofOneWindow),cex=1,srt=45)
				panel.text(nrow(M)/2+5,nrow(M)/2-5,type,cex=1,srt=45)
				panel.text(nrow(M)/2+7,nrow(M)/2-7,CELL,cex=1,srt=45)

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
#------------------------------------------------------------------
OutDir = "tmp"
#------------------------------------------------------------------

MIN_M = 0
MAX_M = 1
# myPalette <- brewer.pal(11, "Spectral")
myPalette <- brewer.pal(9,"YlGnBu")
# FN = paste(OutDir,"/","w_",NumNodesofOneWindow,"_mean.pdf",sep="")
FN = paste(OutDir,"/",CELL,".p1500_d3300.c840.pij.all.pdf",sep="")
cat(FN,"\n")
DrawTriangleHeatmap(FN,NewPIJMat,MIN_M,MAX_M,myPalette,"All")


#------------------------------------------------------------------
# alpha shape
OverlapFN.Common = paste("tmp/",CELL,".p1500_d3300.c840.common.a1.",Type,".txt",sep="")
OverlapFN.New = paste("tmp/",CELL,".p1500_d3300.c840.new.a1.",Type,".txt",sep="")
Overlap.Con.Common = read.table(OverlapFN.Common)
Overlap.Con.New = read.table(OverlapFN.New)
Overlap.Con.All = rbind(Overlap.Con.Common, Overlap.Con.New)
#------------------------------------------------------------------
GetSignificantMat <- function(mat, sig){
	# mat = newmat
	# sig = Overlap.Con.All
	for (i in 1:(NumNodes-2)){
		for (j in (i+2):NumNodes){
			if (length( which(sig[,1]==i & sig[,2]==j )   )==0){
				mat[i,j] = NA
			}
		}
	}
	return(mat)
}
# dev.new()
(FN = paste(OutDir,"/",CELL,".p1500_d3300.c840.",Type,".pij.allsig.pdf",sep=""))
sig.mat.all = GetSignificantMat(NewPIJMat,Overlap.Con.All)
DrawTriangleHeatmap(FN,sig.mat.all,MIN_M,MAX_M,myPalette,"All significant")
# dev.new()
(FN = paste(OutDir,"/",CELL,".p1500_d3300.c840.",Type,".pij.newsig.pdf",sep=""))
sig.mat.new = GetSignificantMat(NewPIJMat,Overlap.Con.New)
DrawTriangleHeatmap(FN,sig.mat.new,MIN_M,MAX_M,myPalette,"New significant")
# dev.new()
(FN = paste(OutDir,"/",CELL,".p1500_d3300.c840.",Type,".pij.commonsig.pdf",sep=""))
sig.mat.common = GetSignificantMat(NewPIJMat,Overlap.Con.Common)
DrawTriangleHeatmap(FN,sig.mat.common,MIN_M,MAX_M,myPalette,"Common significant")
#------------------------------------------------------------------
StripedDiagMatrix <- function(mat,away_s,away_e){
	# mat is a upper triangle matrix
	# mat = sig.mat.all
	# away_s = 10
	# away_e = 54
	NumElement = nrow(mat)
	set_full = seq(1,NumElement)
	# print(away_s)

	for (i in 1:NumElement){
		set_need = intersect(set_full, i + seq(away_s,away_e))
		set_ignore = setdiff(set_full, set_need)
		mat[i,set_ignore] = NA
	}
	tmat = t(mat)
	mat[which(lower.tri(mat))] = tmat[which(lower.tri(tmat))]
	return(mat)
}
#------------------------------------------------------------------
graphics.off()
#------------------------------------------------------------------
away_s = 1
away_e = NumNodes
sig.mat.sel_all = StripedDiagMatrix(sig.mat.all,away_s,away_e)
(FN <- paste(OutDir,"/",CELL,".p1500_d3300.c840.",Type,".pij.all.",away_s,"-",away_e,".pdf",sep=""))
DrawTriangleHeatmap(FN,sig.mat.sel_all,MIN_M,MAX_M,myPalette,paste("Range", away_s,"-",away_e, "significant"))
CI_all = rowSums(sig.mat.sel_all,na.rm=T)
#------------------------------------------------------------------
# away_s = 1
# away_e = 2
# sig.mat.sel_close = StripedDiagMatrix(sig.mat.all,away_s,away_e)
# (FN <- paste(OutDir,"/",CELL,".p1500_d3300.c840.",Type,".pij.all.",away_s,"-",away_e,".pdf",sep=""))
# DrawTriangleHeatmap(FN,sig.mat.sel_close,MIN_M,MAX_M,myPalette,paste("Range", away_s,"-",away_e, "significant"))
# CI_close = rowSums(sig.mat.sel_close,na.rm=T)
# #------------------------------------------------------------------
# away_s = 3
# away_e = NumNodes
# sig.mat.sel_far = StripedDiagMatrix(sig.mat.all,away_s,away_e)
# (FN <- paste(OutDir,"/",CELL,".p1500_d3300.c840.",Type,".pij.all.",away_s,"-",away_e,".pdf",sep=""))
# DrawTriangleHeatmap(FN,sig.mat.sel_far,MIN_M,MAX_M,myPalette,paste("Range", away_s,"-",away_e, "significant"))
# CI_far = rowSums(sig.mat.sel_far,na.rm=T)
#------------------------------------------------------------------
df <- data.frame(
	NodeIndex=1:NumNodes,
	ContactIndex_all = rowSums(sig.mat.sel_all,na.rm=T),
	ContactIndex_close = rowSums(sig.mat.sel_close,na.rm=T),
	ContactIndex_far = rowSums(sig.mat.sel_far,na.rm=T))
df_long <- melt(df,id="NodeIndex")

#------------------------------------------------------------------
p <- ggplot(data=df_long,aes(x=NodeIndex,y=value,color=variable)) +
	theme_bw() +
	geom_line() + geom_point() +
	theme(legend.position = "none") +
	labs(y = "Contact Index") +
	scale_colour_brewer(palette="Set1") +
	scale_x_continuous(breaks=seq(5,NumNodes,by=5))

(FN <- paste(OutDir,"/",CELL,".p1500_d3300.c840.",Type,".pij.all.contactindex.vertical.pdf",sep=""))
pdf(file=FN,width=3, height=5)
p + facet_grid(variable ~.) 
dev.off()

(FN <- paste(OutDir,"/",CELL,".p1500_d3300.c840.",Type,".pij.all.contactindex.horizon.pdf",sep=""))
pdf(file=FN,width=10, height=2)
p + facet_grid(. ~ variable) 
dev.off()

#------------------------------------------------------------------
(FN <- paste(OutDir,"/",CELL,".p1500_d3300.c840.",Type,".pij.all.contactindex.txt",sep=""))
write.table( file=FN,
	paste("#",paste(colnames(df),collapse="\t")), 
	quote=F,row.names=F,col.names=F)
write.table( file =FN,
	apply(df,1,function(x){sprintf("%g\t%.3f\t%.3f\t%.3f",x[1],x[2],x[3],x[4])}),
	quote=F,row.names=F,col.names=F,append=T)
#------------------------------------------------------------------


write.triplet.matrix <- function(m, file="", append=F, fileEncoding=""){
	if (file == ""){
		file <- stdout()
	}
	write.table(paste("#",nrow(m),ncol(m)),file,append=F,quote=F,row.names=F,col.names=F)
	for (i in 1:(nrow(m)-1)){
		for (j in (i+1):ncol(m)){
			if(!is.na(m[i,j]) ){
				if(m[i,j]>0){
					write.table(sprintf("%d\t%d\t%.3e",i,j,m[i,j]),file,append=T,quote=F,row.names=F,col.names=F)
				}
			}
		}
	}
}
#------------------------------------------------------------------
(FN <- paste(OutDir,"/",CELL,".p1500_d3300.c840.",Type,".pij.all.triplet.txt",sep=""))
write.triplet.matrix(file=FN,sig.mat.all)

read.triplet.matrix <- function(file){
	# file = FN
	HeadInfo = read.table(file,comment.char="",nrows=1)
	NumRow = HeadInfo[1,2]
	NumCol = HeadInfo[1,3]

	m = matrix(0,nrow=NumRow,ncol=NumCol)
	DF = read.table(file)
	for (i in 1:nrow(DF)){
		m[DF[i,1],DF[i,2]] = DF[i,3]
		m[DF[i,2],DF[i,1]] = DF[i,3]
	}
	# m = matrix(100,nrow=NumRow,ncol=NumCol)
	# DF = read.table(file)
	# for (i in 1:nrow(DF)){
	# 	if (DF[i,3] != 0){
	# 		cat(i,"\n")
	# 		m[DF[i,2],DF[i,1]] = 1/DF[i,3]
	# 	} else {
	# 		m[DF[i,2],DF[i,1]] = 100
	# 	}
	# }

	return(m)
}
#------------------------------------------------------------------
graphics.off()
dbscan.greater <- function (data, eps, MinPts = 5, scale = FALSE, method = c("hybrid", 
    "raw", "dist"), seeds = TRUE, showplot = FALSE, countmode = NULL) 
{
    distcomb <- function(x, data) {
        data <- t(data)
        temp <- apply(x, 1, function(x) {
            sqrt(colSums((data - x)^2))
        })
        if (is.null(dim(temp))) 
            matrix(temp, nrow(x), ncol(data))
        else t(temp)
    }
    method <- match.arg(method)
    data <- as.matrix(data)
    n <- nrow(data)
    if (scale) 
        data <- scale(data)
    classn <- cv <- integer(n)
    isseed <- logical(n)
    cn <- integer(1)
    for (i in 1:n) {
        if (i %in% countmode) 
            cat("Processing point ", i, " of ", n, ".\n")
        unclass <- (1:n)[cv < 1]
        if (cv[i] == 0) {
            if (method == "dist") {
                reachables <- unclass[data[i, unclass] >= eps]
            }
            else {
                reachables <- unclass[as.vector(distcomb(data[i, 
                  , drop = FALSE], data[unclass, , drop = FALSE])) >= 
                  eps]
            }
            if (length(reachables) + classn[i] < MinPts) 
                cv[i] <- (-1)
            else {
                cn <- cn + 1
                cv[i] <- cn
                isseed[i] <- TRUE
                reachables <- setdiff(reachables, i)
                unclass <- setdiff(unclass, i)
                classn[reachables] <- classn[reachables] + 1
                while (length(reachables)) {
                  if (showplot) 
                    plot(data, col = 1 + cv, pch = 1 + isseed)
                  cv[reachables] <- cn
                  ap <- reachables
                  reachables <- integer()
                  if (method == "hybrid") {
                    tempdist <- distcomb(data[ap, , drop = FALSE], 
                      data[unclass, , drop = FALSE])
                    frozen.unclass <- unclass
                  }
                  for (i2 in seq(along = ap)) {
                    j <- ap[i2]
                    if (showplot > 1) 
                      plot(data, col = 1 + cv, pch = 1 + isseed)
                    if (method == "dist") {
                      jreachables <- unclass[data[j, unclass] >= 
                        eps]
                    }
                    else if (method == "hybrid") {
                      jreachables <- unclass[tempdist[i2, match(unclass, 
                        frozen.unclass)] >= eps]
                    }
                    else {
                      jreachables <- unclass[as.vector(distcomb(data[j, 
                        , drop = FALSE], data[unclass, , drop = FALSE])) >= 
                        eps]
                    }
                    if (length(jreachables) + classn[j] >= MinPts) {
                      isseed[j] <- TRUE
                      cv[jreachables[cv[jreachables] < 0]] <- cn
                      reachables <- union(reachables, jreachables[cv[jreachables] == 
                        0])
                    }
                    classn[jreachables] <- classn[jreachables] + 
                      1
                    unclass <- setdiff(unclass, j)
                  }
                }
            }
        }
        if (!length(unclass)) 
            break
    }
    rm(classn)
    if (any(cv == (-1))) {
        cv[cv == (-1)] <- 0
    }
    if (showplot) 
        plot(data, col = 1 + cv, pch = 1 + isseed)
    out <- list(cluster = cv, eps = eps, MinPts = MinPts)
    if (seeds && cn > 0) {
        out$isseed <- isseed
    }
    class(out) <- "dbscan"
    out
}

get.cluster.dbscan <- function(file,eps,MinPts,cell){
	Node.Dist <- as.dist(read.triplet.matrix(file))
	Node.Mat <- as.matrix(Node.Dist)
	Node.CI <- colSums(Node.Mat)
	dbs <- dbscan.greater(Node.Dist,eps=eps,MinPts=MinPts,method="dist")
	
	df <- data.frame("NodeIndex"=1:length(Node.CI),"CI"=Node.CI, "eps" = eps, "MinPts" = MinPts, "cluster" = as.factor(dbs$cluster), "isseed"=as.factor(if(is.null(dbs$isseed)){F}else{dbs$isseed}),"cell"=as.factor(cell))
	return(df)
}

# CELL = "GM12878"
# CELL = "K562"

# EPS.Vec = c(0.4,0.5,0.6,0.7)
# MinPts.Vec = c(4,5,6)
# EPS.Vec = c(0.6)
# MinPts.Vec = c(4)
#------------------------------------------------------------------
CELL.list = c("GM12878","K562")
CELL.list = c("GM12878")
# CELL.list = c("K562")

EPS.Vec = c(0.4,0.5,0.6,0.7)
MinPts.Vec = c(2,3,4,5,6,7,8,9)

df <- data.frame()
for(k in 1:length(CELL.list)){
	CELL = CELL.list[k]
	(FN <- paste(OutDir,"/",CELL,".p1500_d3300.c840.",Type,".pij.all.triplet.txt",sep=""))
	for (i in 1:length(EPS.Vec)){
		EPS = EPS.Vec[i]
		for (j in 1:length(MinPts.Vec)){
			MINPTS = MinPts.Vec[j]
			df = rbind(df,get.cluster.dbscan(file=FN,eps=EPS,MinPts=MINPTS,cell=CELL))
		}
	}
}

# > display.size <- system("xdpyinfo | grep dimensions", intern = TRUE)
# > display.dpi <- system("xdpyinfo | grep resolution", intern = TRUE)

p <- ggplot(data=df,aes(x=NodeIndex,y=CI,color=cell)) 
p <- p +
	theme_bw() +
	geom_line(size=0.5,colour="gray")  +  
	geom_point(aes(colour = cluster, shape = isseed),size=ifelse(length(EPS.Vec)*length(MinPts.Vec)>6,2,3) ) +
	labs(title=CELL, y = "Contact Index") +
	scale_colour_manual(name="Cluster",values=c("black",brewer.pal(9,"Set1")	)) +
	scale_shape_manual(name="Node",labels=c("border/\nnoise","core"),values=c(1,17)) +
	scale_x_continuous(breaks=seq(5,NumNodes,by=5)) +
	scale_y_continuous(breaks=seq(5,16,by=5))	+
	facet_grid(eps ~ MinPts,labeller = label_both)

dev.new(width=12,height=7)
print(p)
# dev.off()

#------------------------------------------------------------------
# for hinge
CELL.list = c("GM12878","K562")
EPS.Vec = c(0.1)
MinPts.Vec = c(1)

df <- data.frame()
for(k in 1:length(CELL.list)){
	CELL = CELL.list[k]
	(FN <- paste(OutDir,"/",CELL,".p1500_d3300.c840.",Type,".pij.all.triplet.txt",sep=""))
	for (i in 1:length(EPS.Vec)){
		EPS = EPS.Vec[i]
		for (j in 1:length(MinPts.Vec)){
			MINPTS = MinPts.Vec[j]
			df = rbind(df,get.cluster.dbscan(file=FN,eps=EPS,MinPts=MINPTS,cell=CELL))
		}
	}
}

p <- ggplot(data=df,aes(x=NodeIndex,y=CI)) 
p <- p +
	theme_bw() +
	geom_line(size=0.5)  +
	geom_point(aes(colour = cluster, shape = df$isseed),size=ifelse(length(EPS.Vec)*length(MinPts.Vec)>6,2,3) ) +
	labs(title="Hinge", y = "Contact Index") +
	scale_colour_manual(name="Cluster",values=c("black",brewer.pal(9,"Set1")	)) +
	scale_shape_manual(name="Node",labels=c("border/\nnoise","core"),values=c(1,17)) +
	scale_x_continuous(breaks=seq(5,NumNodes,by=5)) +
	scale_y_continuous(breaks=seq(5,16,by=5))	+
	facet_grid(cell ~ eps +MinPts, labeller=label_both)
#------------------------------------------------------------------
library(fpc)
(FN <- paste(OutDir,"/","GM12878",".p1500_d3300.c840.",Type,".pij.all.triplet.txt",sep=""))
GM.Node.Dist <- as.dist(read.triplet.matrix(FN))
(FN <- paste(OutDir,"/","K562",".p1500_d3300.c840.",Type,".pij.all.triplet.txt",sep=""))
K.Node.Dist <- as.dist(read.triplet.matrix(FN))
Diff.Node.Dist <- K.Node.Dist - GM.Node.Dist
Diff.Node.Mat <- as.matrix(Diff.Node.Dist)
Diff.Node.CI <- colSums(Diff.Node.Mat)
EPS.Vec = c(0.1,0.2,0.3)
MinPts.Vec = c(7,8,9)
df <- data.frame()
for (i in 1:length(EPS.Vec)){
	EPS = EPS.Vec[i]
	for (j in 1:length(MinPts.Vec)){
		MINPTS = MinPts.Vec[j]
		Diff.dbs <- dbscan.greater(Diff.Node.Dist,eps=EPS,MinPts=MINPTS,method="dist")
		df = rbind(df,data.frame("NodeIndex"=1:length(Diff.Node.CI),"CI"=Diff.Node.CI, "eps" = EPS, "MinPts" = MINPTS, "cluster" = as.factor(Diff.dbs$cluster), "isseed"=as.factor(if(is.null(Diff.dbs$isseed)){F}else{Diff.dbs$isseed}), "type" = "K-GM" ))
		Diff.dbs <- dbscan(Diff.Node.Dist,eps=-EPS,MinPts=MINPTS,method="dist")
		df = rbind(df,data.frame("NodeIndex"=1:length(Diff.Node.CI),"CI"=Diff.Node.CI, "eps" = -EPS, "MinPts" = MINPTS, "cluster" = as.factor(Diff.dbs$cluster), "isseed"=as.factor(if(is.null(Diff.dbs$isseed)){F}else{Diff.dbs$isseed}), "type" = "GM-K" ))
	}
}

df.subset1 = droplevels(subset(df, type=="K-GM"))
p1 <- ggplot(data=df.subset1,aes(x=NodeIndex,y=CI)) 
p1 <- p1 +
	theme_bw() +
	geom_line(size=0.5,colour="gray")  +  
	geom_point(aes(colour = cluster, shape = isseed),size=ifelse(length(EPS.Vec)*length(MinPts.Vec)>6,2,3) ) +
	labs(title="K-GM", y = "Contact Index Difference ", x= "") +
	scale_colour_manual(name="Cluster",values=c("black",brewer.pal(9,"Set1")	)) +
	scale_shape_manual(name="Node",labels=c("border/\nnoise","core"),values=c(1,17)) +
	scale_x_continuous(breaks=seq(5,NumNodes,by=5)) +
	scale_y_continuous(breaks=seq(5,16,by=5))	+
	facet_grid(eps~ MinPts)

df.subset2 = droplevels(subset(df, type=="GM-K"))
p2 <- ggplot(data=df.subset2,aes(x=NodeIndex,y=CI)) 
p2 <- p2 +
	theme_bw() +
	geom_line(size=0.5,colour="gray")  +  
	geom_point(aes(colour = cluster, shape = isseed),size=ifelse(length(EPS.Vec)*length(MinPts.Vec)>6,2,3) ) +
	labs(title="GM-K", y = "Contact Index Difference ", x="") +
	scale_colour_manual(name="Cluster",values=c("black",brewer.pal(9,"Set1")	)) +
	scale_shape_manual(name="Node",labels=c("border/\nnoise","core"),values=c(1,17)) +
	scale_x_continuous(breaks=seq(5,NumNodes,by=5)) +
	scale_y_continuous(breaks=seq(5,16,by=5))	+
	facet_grid(eps~ MinPts)

dev.new(width=12,height=7)
grid.arrange(p1,p2)

#------------------------------------------------------------------
K.Node.Mat = as.matrix(K.Node.Dist)
K.Node.CI = colSums(K.Node.Mat)
GM.Node.Mat = as.matrix(GM.Node.Dist)
GM.Node.CI = colSums(GM.Node.Mat)
df <- data.frame("NodeIndex"=1:length(Diff.Node.CI),"CI"=GM.Node.CI,"cell"=as.factor("GM"))
df <- rbind(df, data.frame("NodeIndex"=1:length(Diff.Node.CI),"CI"=K.Node.CI,"cell"=as.factor("K")))
df.diff <- data.frame("NodeIndex"=1:length(Diff.Node.CI),"CI"=Diff.Node.CI,"cell"=as.factor("K-GM"))
p <- ggplot(data=df,aes(x=NodeIndex,y=CI,fill=cell,width=.75)) + 
	theme_bw() +
	labs(title="K-GM", y = "Contact Index and Difference ", x="Node Index") +
	geom_line(data=df.diff,aes(x=NodeIndex,y=CI),colour=brewer.pal(9,"Set1")[3],size=1.5) +
	geom_point(data=df.diff,aes(x=NodeIndex,y=CI),size=3,shape=1) +
	geom_bar(stat="identity", position=position_dodge(),alpha=0.5) +
	scale_fill_manual(values=c(brewer.pal(9,"Set1"))) +
	scale_x_continuous(breaks=seq(5,NumNodes,by=5)) 
dev.new(width=12,height=7)
p	
#------------------------------------------------------------------
#------------------------------------------------------------------
#------------------------------------------------------------------
