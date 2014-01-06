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
  #CELL = "GM12878"
}else{
  CELL = args[1]
}

#------------------------------------------------------------------
Type="con"
method = "cut.all"
#------------------------------------------------------------------
ConDir=paste("/Users/yunxu/workspace/projects/chromatin/result/analysis/ENm008/chain/",method,"/par.1_1_1_1/",CELL,".p1500_d3300",sep="")
# ConDir=paste("/dump/yunxu/working/projects.new/chromatin/result/analysis/ENm008/chain/",method,"/par.1_1_1_1/",CELL,".p1500_d3300",sep="")
#------------------------------------------------------------------

NumNodes = 54
DistCutoff = 840
WeightDistPath = paste("/Users/yunxu/workspace/projects/chromatin/result/analysis/ENm008/chain/weight.dist/par.1_1_1_1/",CELL,".p1500_d3300",sep="")


#------------------------------------------------------------------
Pair_Vec = t(combn(1:NumNodes,2))
Pair_Vec = Pair_Vec[Pair_Vec[,2]-Pair_Vec[,1]!=1,]

#------------------------------------------------------------------
# Get Pij of chain
FileList = apply(Pair_Vec,1,function(x){
        paste(WeightDistPath,"/",paste(x,collapse="_"), ".txt",sep="")})

SampleFile = FileList[1]
SampleData = read.table(SampleFile)
NRow = nrow(SampleData)
NCol = length(FileList)
Weight_Vec = SampleData[,2]

DistMat = matrix(0, nrow=NRow, ncol=NCol)
SignMat = DistMat
for (i in 1:NCol){
  FN = FileList[i]
  WeightDistData = read.table(FN)
  DistMat[,i] = WeightDistData[,1]
}
SignMat[DistMat<DistCutoff] = 1
Pij = apply(SignMat,2,function(x){
  SignInd = which(x==1)        
  sum(Weight_Vec[SignInd])/sum(Weight_Vec)
        })

PijMat = matrix(NA,nrow=NumNodes,ncol=NumNodes)
for (i in 1:(NumNodes-2)){
  for (j in (i+2):NumNodes){
    idx = which(Pair_Vec[,1] == i & Pair_Vec[,2] == j)
    PijMat[i,j] = Pij[idx]
  }
}

# PPij = apply(SignMat,2,function(x){
#   SignInd = which(x==1)        
#   sum(Weight_Vec[SignInd])
#         })
# PPijMat = matrix(0,nrow=NumNodes, ncol=NumNodes)
# for (i in 1:(NumNodes-2)){
#   for (j in (i+2):NumNodes){
#     idx = which(Pair_Vec[,1] == i & Pair_Vec[,2] == j)
#     PPijMat[i,j] = PPij[idx]
#   }
# }
# get significant contact
#------------------------------------------------------------------
SignificantPath = "/Users/yunxu/workspace/projects/chromatin/result/analysis/ENm008/chain/sts/pij/"
OverlapFN.Common = paste(SignificantPath,CELL,".p1500_d3300.c840.common.a1.",Type,".txt",sep="")
OverlapFN.New = paste(SignificantPath,CELL,".p1500_d3300.c840.new.a1.",Type,".txt",sep="")
Overlap.Con.Common = read.table(OverlapFN.Common)
Overlap.Con.New = read.table(OverlapFN.New)
Overlap.Con.All = rbind(Overlap.Con.Common, Overlap.Con.New)

GetSignificantMat <- function(mat, sig){
#  mat = PijMat
#  sig = Overlap.Con.All
  for (i in 1:(NumNodes-2)){
    for (j in (i+2):NumNodes){
      # in significant calculation, I use <= may miss the pij = 0
      # so here consider pij=0 situation
      if (length( which(sig[,1]==i & sig[,2]==j ) )==0){
        mat[i,j] = NA
      }
    }
  }
  return(mat)
}
PijMat_Sig = GetSignificantMat(PijMat,Overlap.Con.All)

#ConIndex = apply(t(PijMat_Sig)+PijMat_Sig, 1,function(x){ sum(x,na.rm=T)})
#------------------------------------------------------------------
# Drawing triangle
#------------------------------------------------------------------
#----------------------------------------------------------------------------
All_Points_Pos = c(1,5693,15091,18344,29756,44231,50868,55911,64056,74448,88084,95256,100530,104687,109838,123474,131220,134334,147970,161606,167103,171769,185994,189074,203353,217802,225341,238977,247214,260850,274486,277942,289198,303867,310528,314538,327240,334889,348525,352385,360924,372670,380160,393796,407432,418222,421291,433293,445126,454365,468001,483412,496513,499411)
GetGeneNameAndRange <- function(){
  FN = "data/ENm008/circos/text.genes.txt"
  GeneRangeName = read.table(FN)[2:4]
  PosToName = NULL
  for (i in 1:length(All_Points_Pos)){
    pos = All_Points_Pos[i]
    str = ""
    for (j in 1:nrow(GeneRangeName)){
      if (  GeneRangeName[j,1] <= pos && GeneRangeName[j,2] >= pos ){
        if (str == ""){
          str = as.character(GeneRangeName[j,3])
        } else{
          str = paste(str, as.character(GeneRangeName[j,3]),sep=" / ")
        }
      }
    }
    # cat (i,"\t",str,"\n")
    PosToName = rbind(PosToName,c(i,str))
  }
  return(PosToName)
}

#----------------------------------------------------------------------------
PosToName = GetGeneNameAndRange()
PosToName[! (PosToName[,2] %in% c("HS48 / HS40","HS40",
                                  "HS10","HS10 / HBZ/HS8",
                                  "HBM/HBA2/HBA1",
                                  "HBM/HBA2/HBA1 / HBQ1")),2]= ""
PosToName = PosToName[seq(5,50,by=5),]
xlabel = apply(PosToName,1,function(x){paste(x[2],x[1])})
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
  # remove zero 
  Mat[Mat==0] = NA
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
#        xlabel.fmt = sprintf("_% 23s",xlabel)
        #for (i in 1:nrow(M)){
##          panel.text(i+0.25,i-0.5,xlabel[i],cex=.5,srt=90,adj=c(1,NA))
          #panel.text(i,i+0.5,xlabel[i],cex=.5,srt=90,adj=c(1,NA))
        #}

        for (i in 1:nrow(PosToName)){
          panel.text(
            as.numeric(PosToName[i,1]),
            as.numeric(PosToName[i,1])+0.5,PosToName[i,1],cex=1,srt=90,adj=c(1,NA))
        }
        # panel.text(nrow(M)/2+3,nrow(M)/2-3,paste("Window size =",NumNodesofOneWindow),cex=1,srt=45)
        #panel.text(nrow(M)/2+6,nrow(M)/2-6,type,cex=1,srt=45)
        #panel.text(nrow(M)/2+8,nrow(M)/2-8,CELL,cex=1,srt=45)

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
  legend.vp <- viewport(  x = unit(2.5,"lines"),  y = unit(0.7,"npc"),
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

OutDir = "tmp1"
MIN_M = 0
MAX_M = 1
# myPalette <- brewer.pal(11, "Spectral")
myPalette <- brewer.pal(9,"YlGnBu")
myPalette <- brewer.pal(9,"OrRd")
FN = paste(OutDir,"/",CELL,".p1500_d3300.c840.pij.all.full.pdf",sep="")
cat(FN,"\n")
DrawTriangleHeatmap(FN,PijMat,MIN_M,MAX_M,myPalette,"All Pij")

FN = paste(OutDir,"/",CELL,".p1500_d3300.c840.pij.all.sig.pdf",sep="")
cat(FN,"\n")
DrawTriangleHeatmap(FN,PijMat_Sig,MIN_M,MAX_M,myPalette,"All Significant")

PijMat_Sig_New = GetSignificantMat(PijMat_Sig,Overlap.Con.New)
FN = paste(OutDir,"/",CELL,".p1500_d3300.c840.pij.new.sig.pdf",sep="")
cat(FN,"\n")
DrawTriangleHeatmap(FN,PijMat_Sig_New,MIN_M,MAX_M,myPalette,"New Significant")

PijMat_Sig_Common = GetSignificantMat(PijMat_Sig,Overlap.Con.Common)
FN = paste(OutDir,"/",CELL,".p1500_d3300.c840.pij.common.sig.pdf",sep="")
cat(FN,"\n")
DrawTriangleHeatmap(FN,PijMat_Sig_Common,MIN_M,MAX_M,myPalette,"Overlap Significant")


#------------------------------------------------------------------
Get_PijConnection_CircosStr <- function(FileName,Mat){
  index.vec = which(!is.na(Mat) & Mat != 0,arr.ind=T)

  idx = 1
  str_circos = apply(index.vec,1,function(x){
    str_head = sprintf("segdup%04d hs16",idx);
    str1 = paste(str_head,
      All_Points_Pos[x[1]],All_Points_Pos[x[1]],
      sprintf("pij=%.2f",Mat[x[1],x[2]]))
    str2 = paste(str_head,
      All_Points_Pos[x[2]],All_Points_Pos[x[2]],
      sprintf("pij=%.2f",Mat[x[1],x[2]]))

    idx <<- idx + 1;
    paste(str1,str2,sep="\n")
    }
  )
  cat(CELL,"\n")
  write.table(file=FileName, 
    paste(str_circos,collapse="\n"),col.names=F,row.names=F,quote=F)
}

(FN = paste("/tmp/",CELL,".Pred.All.Sig.Pij.Circos.txt",sep=""))
Get_PijConnection_CircosStr(FN,PijMat_Sig)
