#!/usr/bin/R
rm(list=ls())
WorkingDir="/Users/yunxu/workspace/projects/chromatin/result/analysis/ENm008/chain/par.1_1_1_1/"
WorkingDir="/dump/yunxu/working/projects.new/chromatin/result/analysis/ENm008/chain/par.1_1_1_1/"
SampleSize=c("GM12878.p1500_d3300","K562.p1500_d3300")

NodeNodeFile = "nodenodelist.txt"
NodeNodeVec = read.table(NodeNodeFile)[,1:2]

WeightDistGrp = NULL
ConVec = NodeNodeVec[which(NodeNodeVec[,1]==20 | NodeNodeVec[,2]==20),]
# ConVec = ConVec[1:2,]

for (i in 1:nrow(ConVec)){
  cat(i,"/",nrow(ConVec),"\n")
  flush.console()

  FN = paste(WorkingDir,"/",SampleSize[1],"/",ConVec[i,1],"_",ConVec[i,2],".txt",sep="")
  WD = read.table(FN)

  # WD = WD[1:100,]
  Dist = NULL
  for (j in 1:nrow(WD)){
    Dist = c(Dist, rep(WD[j,1],WD[j,2]))
  }
  WeightDistGrp = rbind(WeightDistGrp, 
    data.frame(
       dist=Dist,
       cell=rep(factor("GM"),length(Dist)),
       con=rep(factor(paste(ConVec[i,1],"_",ConVec[i,2],sep="")),length(Dist)))
    )


  FN = paste(WorkingDir,"/",SampleSize[2],"/",ConVec[i,1],"_",ConVec[i,2],".txt",sep="")
  WD = read.table(FN)

  # WD = WD[1:100,]
  Dist = NULL
  for (j in 1:nrow(WD)){
    Dist = c(Dist, rep(WD[j,1],WD[j,2]))
  }
  WeightDistGrp = rbind(WeightDistGrp, 
    data.frame(
       dist=Dist,
       cell=rep(factor("K"),length(Dist)),
       con=rep(factor(paste(ConVec[i,1],"_",ConVec[i,2],sep="")),length(Dist)))
    )
}
# convert A to nm
WeightDistGrp[,"dist"] = WeightDistGrp[,"dist"] / 10



GM.ConFile = "/tmp/GM12878.p1500_d3300.c840.a5.sis.node.dst"
K.ConFile = "/tmp/K562.p1500_d3300.c840.a5.sis.node.dst"
GM.ConVec = read.table(GM.ConFile)
K.ConVec = read.table(K.ConFile)

GM.ConVec = GM.ConVec[which(GM.ConVec[,1]==20 | GM.ConVec[,2]==20),]
K.ConVec = K.ConVec[which(K.ConVec[,1]==20 | K.ConVec[,2]==20),]

GM.factor = factor(apply(GM.ConVec,1,function(x){paste(x[1],x[2],sep="_")}))
K.factor = factor(apply(K.ConVec,1,function(x){paste(x[1],x[2],sep="_")}))
whole.factor = levels(WeightDistGrp[,3])

dev.new(width=12, height=6)
xlabels = levels(WeightDistGrp[,3])
boxplot(dist ~ con, data = WeightDistGrp,
        boxwex = 0.25,  at = 1:length(xlabels) - 0.2,
        subset = cell == "GM", col = "yellow",
        main = "Alpha gene node pair distance boxplot",
        xlab = "Alpha gene node pair",
        ylab = "Distance (nm)",
        outline=F,
        xlim = c(1,length(xlabels)),
        ylim = range(WeightDistGrp[,1]),
        axes=F)
boxplot(dist ~ con, data = WeightDistGrp, add=TRUE,
        boxwex = 0.25, at = 1:length(levels(WeightDistGrp[,3])) + 0.2,
        subset = cell == "K", col = "orange",
        outline=F,
        axes=F)
box()
axis(1,at=1:length(xlabels),labels=FALSE)
axis(2)
text(1:length(xlabels),par("usr")[3]-50,labels = xlabels, srt = 45, pos = 2, xpd = TRUE,offset=c(-0.2,1))
legend("topleft",c("GM12878","K562"),fill=c("yellow","orange"))

dev.new(width=12, height=6)
plot(0,xlim=c(1,length(xlabels)),ylim=range(WeightDistGrp[,1]),type="n",
     xlab="",ylab="",axes=F)
boxplot(dist ~ con, data = WeightDistGrp,
        boxwex = 0.25,  at = 1:length(xlabels) - 0.2, add=TRUE,
        subset = (cell == "GM" & con %in% GM.factor), col = "yellow",
        main = "Alpha gene node pair distance boxplot",
        xlab = "Alpha gene node pair",
        ylab = "Distance (nm)",
        outline=F,
        xlim = c(1,length(xlabels)),
        ylim = range(WeightDistGrp[,1]),
        axes=F)
boxplot(dist ~ con, data = WeightDistGrp, add=TRUE,
        boxwex = 0.25, at = 1:length(levels(WeightDistGrp[,3])) + 0.2,
        subset = (cell == "K" & con %in% K.factor), col = "orange",
        outline=F,
        axes=F)
box()
axis(1,at=1:length(xlabels),labels=FALSE)
axis(2)
text(1:length(xlabels),par("usr")[3]-50,labels = xlabels, srt = 45, pos = 2, xpd = TRUE,offset=c(-0.2,1))
legend("topleft",c("GM12878","K562"),fill=c("yellow","orange"))
