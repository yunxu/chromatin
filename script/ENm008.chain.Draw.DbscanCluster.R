#!/usr/bin/R
rm(list=ls())
#------------------------------------------------------------------
library(grid)
library(gridExtra)
library(lattice)
library(RColorBrewer)
library(reshape)
library(ggplot2)
library(fpc)
#------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (Sys.info()["sysname"] == "Darwin"){
  CELL = "K562"
  CELL = "GM12878"
}else{
  CELL = args[1]
}

#------------------------------------------------------------------
NumNodes = 54
DistCutoff = 840
Type="con"
method = "cut.all"
#------------------------------------------------------------------
ConDir=paste("/Users/yunxu/workspace/projects/chromatin/result/analysis/ENm008/chain/",method,"/par.1_1_1_1/",CELL,".p1500_d3300",sep="")
# ConDir=paste("/dump/yunxu/working/projects.new/chromatin/result/analysis/ENm008/chain/",method,"/par.1_1_1_1/",CELL,".p1500_d3300",sep="")
#----------------------------------------------------------------------------
GetGeneNameAndRange <- function(){
  All_Points_Pos = c(1,5693,15091,18344,29756,44231,50868,55911,64056,74448,88084,95256,100530,104687,109838,123474,131220,134334,147970,161606,167103,171769,185994,189074,203353,217802,225341,238977,247214,260850,274486,277942,289198,303867,310528,314538,327240,334889,348525,352385,360924,372670,380160,393796,407432,418222,421291,433293,445126,454365,468001,483412,496513,499411)
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
xlabel = apply(PosToName,1,function(x){paste(x[2],x[1])})
xlabel = seq(1,54,by=1)
xlabel = seq(5,54,by=5)

#------------------------------------------------------------------
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
#------------------------------------------------------------------
Get_PijMat_ByCell <- function(cell){

  WeightDistPath = paste("/Users/yunxu/workspace/projects/chromatin/result/analysis/ENm008/chain/weight.dist/par.1_1_1_1/",cell,".p1500_d3300",sep="")
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
  
  #------------------------------------------------------------------
  SignificantPath = "/Users/yunxu/workspace/projects/chromatin/result/analysis/ENm008/chain/sts/pij/"
  OverlapFN.Common = paste(SignificantPath,cell,".p1500_d3300.c840.common.a1.",Type,".txt",sep="")
  OverlapFN.New = paste(SignificantPath,cell,".p1500_d3300.c840.new.a1.",Type,".txt",sep="")
  Overlap.Con.Common = read.table(OverlapFN.Common)
  Overlap.Con.New = read.table(OverlapFN.New)
  Overlap.Con.All = rbind(Overlap.Con.Common, Overlap.Con.New)
  
  PijMat_Sig = GetSignificantMat(PijMat,Overlap.Con.All)
  PijMat_Sig[PijMat_Sig==0] = NA
  return(PijMat_Sig)
}

PijMat_Sig_GM = Get_PijMat_ByCell("GM12878")
PijMat_Sig_K = Get_PijMat_ByCell("K562")



#------------------------------------------------------------------
# Drawing triangle
#------------------------------------------------------------------
#------------------------------------------------------------------
OutDir = "tmp1"
#------------------------------------------------------------------
#------------------------------------------------------------------
#------------------------------------------------------------------
#------------------------------------------------------------------
#------------------------------------------------------------------
#------------------------------------------------------------------
#------------------------------------------------------------------
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

#------------------------------------------------------------------
get.cluster.dbscan <- function(pijmat,eps,MinPts,cell){
  # pijmat is upper triangle matrix
  #pijmat <- PijMat_Sig
  pijmat[is.na(pijmat)] = 0
  Node.Mat <- pijmat + t(pijmat)
  Node.Dist <- as.dist(Node.Mat)
  Node.CI <- colSums(Node.Mat)
  dbs <- dbscan.greater(Node.Dist,eps=eps,MinPts=MinPts,method="dist")
  
  df <- data.frame("NodeIndex"=1:length(Node.CI),"CI"=Node.CI, "eps" = eps, "MinPts" = MinPts, "cluster" = as.numeric(dbs$cluster), "isseed"=as.factor(if(is.null(dbs$isseed)){F}else{dbs$isseed}),"cell"=as.factor(cell))
  return(df)
}

#------------------------------------------------------------------

MyPal = brewer.pal(9,"Set1")
EPS = 0.6
MINPTS = 7
df_GM = get.cluster.dbscan(PijMat_Sig_GM, eps = EPS, MinPts=MINPTS, cell="GM12878")
df_K = get.cluster.dbscan(PijMat_Sig_K, eps = EPS, MinPts=MINPTS, cell="K562") 
df_K[df_K["cluster"]==1,"cluster"] =2


df = rbind(df_GM, df_K)
df$cluster <- as.factor(df$cluster)

gp <- ggplot(data=df, aes(x=NodeIndex,y=CI,color=cluster))
gp = gp +
  theme_bw() +
  theme(
    text = element_text(size=25),
    axis.title.x=element_text(vjust=-1),
    #axis.text.x=element_text(angle=30, hjust=1,vjust=1)
    legend.position = "none"
  ) +
  #scale_x_discrete(label=xlabel) +
  scale_x_discrete(breaks=xlabel,label=xlabel) +
  scale_y_discrete(breaks=seq(0,15,by=5),label=seq(0,15,by=5)) +
  coord_cartesian(ylim=c(-1,17)) +
  labs(x= paste("Primer sites"),y = "Contact Index")

gp = gp + geom_line(size=1, colour="grey") +
  #geom_point(aes(colour = as.factor(cluster), shape = isseed),size=2 ) +
  #geom_point(data=subset(df,cluster %in% c(1,2)),aes(colour = cluster, shape = isseed),size=2.5 ) +
  geom_point(aes(colour = as.factor(cluster), shape = isseed),size=5 ) +
  geom_point(data=subset(df,cluster %in% c(1,2)),aes(colour = cluster, shape =
                                                     isseed),size=4.5 ) +
  scale_colour_manual(name="Cluster",values=c("black",MyPal) )+
  scale_shape_manual(name="Primer Site",labels=c("border/\nnoise","core"),values=c(1,17)) 

#dev.new(width=17,height=4)
dev.new(width=14,height=4.5)
  gp + facet_grid(cell~.)
dev.copy2eps(file="/tmp/dbcluster.eps")
#postscript("/tmp/dbcluster.eps",width=14,height=4.5)
dev.off()
