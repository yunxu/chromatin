#!/usr/bin/R

rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)

if (Sys.info()["sysname"] == "Darwin"){
	NodeInd = args[1]
	NodeInd = 21
	WorkingDir="/Users/yunxu/workspace/projects/chromatin/result/analysis/ENm008/chain/weight.dist/par.1_1_1_1/"
	
}else{
	NodeInd = args[1]
	WorkingDir="/dump/yunxu/working/projects.new/chromatin/result/analysis/ENm008/chain/weight.dist/par.1_1_1_1/"	
}
NodeNodeFile = "data/analysis/ENm008/nodenodelist.txt"
GM.ConFile = "data/analysis/ENm008/GM12878.p1500_d3300.c840.a5.sis.node.dst"
K.ConFile = "data/analysis/ENm008/K562.p1500_d3300.c840.a5.sis.node.dst"

# WorkingDir="/Users/yunxu/workspace/projects/chromatin/result/analysis/ENm008/chain/weight.dist/par.1_1_1_1/"
#WorkingDir="/dump/yunxu/working/projects.new/chromatin/result/analysis/ENm008/chain/par.1_1_1_1/"
# WorkingDir="par.1_1_1_1/"
OutDir = "/tmp/test"
if (!file_test("-d",OutDir)){dir.create(OutDir)}

SampleSize=c("GM12878.p1500_d3300","K562.p1500_d3300")
NodeNodeVec = read.table(NodeNodeFile)[,1:2]

GM.ConVec = read.table(GM.ConFile)
K.ConVec = read.table(K.ConFile)

GM.ConVec[,3] = GM.ConVec[,3] / 10
K.ConVec[,3] = K.ConVec[,3] / 10

ConVec = NodeNodeVec[which(NodeNodeVec[,1]==NodeInd | NodeNodeVec[,2]==NodeInd),]
GM_NodeVec = GM.ConVec[which(GM.ConVec[,1]==NodeInd | GM.ConVec[,2]==NodeInd),]
K_NodeVec = K.ConVec[which(K.ConVec[,1]==NodeInd | K.ConVec[,2]==NodeInd),]

GM.Ind = as.numeric(apply(GM_NodeVec,1,function(x){which(ConVec[,1]==x[1] & ConVec[,2]==x[2])}))
K.Ind = as.numeric(apply(K_NodeVec,1,function(x){which(ConVec[,1]==x[1] & ConVec[,2]==x[2])}))
# ConVec = ConVec[1:5,]

# ------------------------------------------------------------------
FN = paste(WorkingDir,"/",SampleSize[1],"/",ConVec[1,1],"_",ConVec[1,2],".txt",sep="")
numRow = sum(read.table(FN)[,2])
WeightDistGrp1 = matrix(NA,nrow=numRow,ncol=nrow(ConVec))
for (i in 1:nrow(ConVec)){
#  cat(i,"/",nrow(ConVec),"\n")
#  flush.console()
  FN = paste(WorkingDir,"/",SampleSize[1],"/",ConVec[i,1],"_",ConVec[i,2],".txt",sep="")
  WD = read.table(FN)

	WeightDistGrp1[,i] = rep(WD[,1],WD[,2])
}
WeightDistGrp1 = WeightDistGrp1/10

# ------------------------------------------------------------------
FN = paste(WorkingDir,"/",SampleSize[2],"/",ConVec[1,1],"_",ConVec[1,2],".txt",sep="")
numRow = sum(read.table(FN)[,2])
WeightDistGrp2 = matrix(NA,nrow=numRow,ncol=nrow(ConVec))
for (i in 1:nrow(ConVec)){
#  cat(i,"/",nrow(ConVec),"\n")
#  flush.console()
  FN = paste(WorkingDir,"/",SampleSize[2],"/",ConVec[i,1],"_",ConVec[i,2],".txt",sep="")
  WD = read.table(FN)

	WeightDistGrp2[,i] = rep(WD[,1],WD[,2])
}
WeightDistGrp2 = WeightDistGrp2/10
# ------------------------------------------------------------------



BoxplotWeightDist_Sel <- function(){
	# dev.new(width=12, height=6)
	# FN = paste(OutDir,"/",NodeInd,"_in.pdf",sep="")
	# pdf(FN,width=12,height=6)
	X_AtPos_1 = GM.Ind - 0.2
	X_AtPos_2 = K.Ind + 0.2
	offset = 20
	
	plot(0,
		xlim=c(1,length(xlabels)),
		# ylim=range(WeightDistGrp1,WeightDistGrp2),
		ylim = c(0,800),
		type="n", xlab="",ylab="",axes=F)
	if (length(GM.Ind)!=0)	{
		boxplot(WeightDistGrp1[,GM.Ind],
			boxwex = 0.25,  
			at = GM.Ind - 0.2, add=TRUE,
			col = "yellow",
			# main = paste("Node ",NodeInd," pair distance boxplot"),
			# xlab = paste("Node ",NodeInd," pair"),
			# ylab = "Distance (nm)",
			outline=F,
			axes=F)
		
		#     MedianValue_1 = sprintf("%1.0f",apply(WeightDistGrp1[,GM.Ind,drop=FALSE],2,median))
		# Y_AtPos_1 = apply(WeightDistGrp1[,GM.Ind,drop=FALSE],2,function(x){
		# 	y = quantile(x,prob=c(.25,.75))
		# 	y = y[2] + 1.5*(y[2]-y[1]) + offset
		# 	})
		# text(X_AtPos_1,Y_AtPos_1,MedianValue_1,cex=.5)
		
    MeanValue_1 = sprintf("%1.0f",apply(WeightDistGrp1[,GM.Ind,drop=FALSE],2,mean))
		Y_AtPos_1 = apply(WeightDistGrp1[,GM.Ind,drop=FALSE],2,function(x){
			y = quantile(x,prob=c(.25,.75))
			y = y[2] + 1.5*(y[2]-y[1]) + 3 *offset
			})
		text(X_AtPos_1,Y_AtPos_1,MeanValue_1,cex=.5,col="red")
		}
	if (length(K.Ind)!=0){
		boxplot(WeightDistGrp2[,K.Ind], add=TRUE,
			boxwex = 0.25, 
			at = K.Ind + 0.2,
			col = "orange",
			# xlab = paste("Node ",NodeInd," pair"),
			# ylab = "Distance (nm)",
			outline=F,
			axes=F)

			# MedianValue_2 = sprintf("%1.0f",apply(WeightDistGrp2[,K.Ind,drop=FALSE],2,median))
			# Y_AtPos_2 = apply(WeightDistGrp2[,K.Ind,drop=FALSE],2,function(x){
			# 	y = quantile(x,prob=c(.25,.75))
			# 	y = y[2] + 1.5*(y[2]-y[1]) + offset
			# 	})
			# text(X_AtPos_2,Y_AtPos_2,MedianValue_2,cex=.5)	
			
      MeanValue_2 = sprintf("%1.0f",apply(WeightDistGrp2[,K.Ind,drop=FALSE],2,mean))
			Y_AtPos_2 = apply(WeightDistGrp2[,K.Ind,drop=FALSE],2,function(x){
				y = quantile(x,prob=c(.25,.75))
				y = y[2] + 1.5*(y[2]-y[1]) + 3 * offset
				})
			text(X_AtPos_2,Y_AtPos_2,MeanValue_2,cex=.5,col="red")	
		}
	box()
	mtext("Distance (nm)",side=2,line=3)
	axis(1,at=1:length(xlabels),labels=FALSE)
	axis(2)
	text(1:length(xlabels),par("usr")[3]-50,labels = xlabels, srt = 45, pos = 2, xpd = TRUE,offset=c(-0.2,1))
	legend("topleft",c("GM12878","K562"),fill=c("yellow","orange"))
}
BoxplotWeightDist_All <- function(){
	# dev.new(width=12, height=6)
	# FN = paste(OutDir,"/",NodeInd,"_all.pdf",sep="")
	# pdf(FN,width=12,height=6)
	X_AtPos_1 = 1:length(xlabels) - 0.2
	X_AtPos_2 = 1:length(xlabels) + 0.2
	offset = 20
	boxplot(WeightDistGrp1,
		boxwex = 0.25,  
		at = X_AtPos_1,
		col = "yellow",
		# main = paste("Node ",NodeInd," pair distance boxplot"),
		xlab = paste("Node ",NodeInd," pair"),
		ylab = "Distance (nm)",
		outline=F,
		cex = 0.5,
		xlim = c(1,length(xlabels)+0.2),
		# ylim = range(WeightDistGrp1,WeightDistGrp2),
		ylim = c(0,800),
		axes=F)
	boxplot(WeightDistGrp2, add=TRUE,
		boxwex = 0.25, 
		at = X_AtPos_2,
		col = "orange",
		outline=F,
		axes=F)	
	box()
	axis(1,at=1:length(xlabels),labels=FALSE)
	axis(2)
	text(1:length(xlabels),par("usr")[3]-50,labels = xlabels, srt = 45, pos = 2, xpd =
	 TRUE,offset=c(-0.2,1))
	legend("topleft",c("GM12878","K562"),fill=c("yellow","orange"))
	
	# MedianValue_1 = sprintf("%1.0f",apply(WeightDistGrp1,2,median))
	# Y_AtPos_1 = apply(WeightDistGrp1,2,function(x){
	# 	y = quantile(x,prob=c(.25,.75))
	# 	y = y[2] + 1.5*(y[2]-y[1]) + offset
	# 	})
	# text(X_AtPos_1,Y_AtPos_1,MedianValue_1,cex=.5)
	
  MeanValue_1 = sprintf("%1.0f",apply(WeightDistGrp1,2,mean))
	Y_AtPos_1 = apply(WeightDistGrp1,2,function(x){
		y = quantile(x,prob=c(.25,.75))
		y = y[2] + 1.5*(y[2]-y[1]) + 3 * offset
		})
	text(X_AtPos_1,Y_AtPos_1,MeanValue_1,cex=.5,col="red")
	
	# MedianValue_2 = sprintf("%1.0f",apply(WeightDistGrp2,2,median))
	# Y_AtPos_2 = apply(WeightDistGrp2,2,function(x){
	# 	y = quantile(x,prob=c(.25,.75))
	# 	y = y[2] + 1.5*(y[2]-y[1]) + offset
	# 	})
	# text(X_AtPos_2,Y_AtPos_2,MedianValue_2,cex=.5)	
	
  MeanValue_2 = sprintf("%1.0f",apply(WeightDistGrp2,2,mean))
	Y_AtPos_2 = apply(WeightDistGrp2,2,function(x){
		y = quantile(x,prob=c(.25,.75))
		y = y[2] + 1.5*(y[2]-y[1]) + 3 * offset
		})
	text(X_AtPos_2,Y_AtPos_2,MeanValue_2,cex=.5,col="red")	
}

OneFigure <- function(){
	FN = paste(OutDir,"/all_",NodeInd,".pdf",sep="")
  cat(FN,"\n")
	pdf(FN,width=12,height=8)
	# pdf(FN)
	
	par(mfrow=c(2,1),mar=c(5.1,4.1,0.5,0.5),oma=c(0,0,2,0))
	BoxplotWeightDist_Sel()
	# BoxplotWeightDist_Sel()
	BoxplotWeightDist_All()
	mtext(paste("Node ",NodeInd," pair distance boxplot"), cex=1.2, side=3, outer=TRUE) 
	dev.off()
}

xlabels = apply(ConVec,1,function(x){paste(x[1],x[2],sep="_")})
OneFigure()
# ------------------------------------------------------------------

FN = paste(OutDir,"/GM_",NodeInd,"_sts.txt",sep="")
write.table(
	file=FN,
	cbind(
		sprintf("%10.3f",GM.ConVec[GM.Ind,3]),
		sprintf("%10.3f",apply(WeightDistGrp1,2,median)[GM.Ind])),
	quote=F,row.names=F,col.names=F,sep="\t")

FN = paste(OutDir,"/K_",NodeInd,"_sts.txt",sep="")
write.table(
	file=FN,
	cbind(
		sprintf("%10.3f",K.ConVec[K.Ind,3]),
		sprintf("%10.3f",apply(WeightDistGrp2,2,median)[K.Ind])),
	quote=F,row.names=F,col.names=F,sep="\t")

FN = paste(OutDir,"/","mean.sd.dist_",NodeInd,".txt",sep="")
write.table(file=FN,"#Node Node GM12878(mean) GM12878(sd) K562(mean) K562(sd)",	quote=F,row.names=F,col.names=F,sep="\t")
write.table(file=FN,
	cbind(ConVec, sprintf("%1.0f",apply(WeightDistGrp1,2,mean)),  sprintf("%1.0f",apply(WeightDistGrp1,2,sd)),  		
		sprintf("%1.0f",apply(WeightDistGrp2,2,mean)),sprintf("%1.0f",apply(WeightDistGrp2,2,sd))),
	quote=F,row.names=F,col.names=F,sep="\t",
	append=T)
