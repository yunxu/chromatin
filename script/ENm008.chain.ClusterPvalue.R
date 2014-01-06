#!/usr/bin/R
#------------------------------------------------------------------
# ENm008.chain.ClusterPvalue.R
#------------------------------------------------------------------
rm(list=ls())
# library(fpc)
#------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
ListFN <- args[1]
ListFN <- "/tmp/subaae"
ListFN <- "tmp/sub/subaae"
Type = "r270"
method = "alpha"
NumNodes = 54
MaxLogWeight = 326.556
OutDir="tmp/cluster/"
if ( ! file.exists(OutDir)){dir.create(OutDir)}
#------------------------------------------------------------------
if (Sys.info()["sysname"] == "Darwin"){
	DATADIR = "/Users/yunxu/workspace/projects/chromatin/"
	ANALYSISDIR = "/Users/yunxu/workspace/projects/chromatin/"
}else{
	DATADIR = "/dump/yunxu/working/projects.new/chromatin/"
	ANALYSISDIR = "/dump/yunxu/working/projects.new/chromatin/"
}
PTSDIR = paste(DATADIR,"result/ENm008/ensemble2/p1500_d3300/",sep="")
OVERLAPDIR = paste(ANALYSISDIR,"result/analysis/ENm008/ensemble2/alpha/",sep="")
#------------------------------------------------------------------

PTSFiles = list.files(path=PTSDIR,pattern=".pts",full.names=T,recursive=T)
PTSFiles = as.character(read.table(ListFN)[,1])
NumFiles = length(PTSFiles)
BaseNames = basename(sub("[.][^.]*$","",PTSFiles))
SubDirs = basename(dirname(PTSFiles))
OverlapFiles = paste(OVERLAPDIR,SubDirs,"/",BaseNames,".",Type,".05.a.overlap",sep="")
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
GetWeightDensityCluster_GrpGrp <- function(ptsfile,confile,eps,MinPts,GrpList){
	# ptsfile = PTSFiles[77]
	# confile = OverlapFiles[77]
	# eps = EPS
	# MinPts = MINPTS
	# Grp1 = 14
	# Grp2 = 20
	# GrpList = list(list(Grp1),list(Grp2))
	LogWeight = read.table(ptsfile,comment.char="",nrow=1)[,3]
	ConMatrix = matrix(0,nrow=NumNodes,ncol=NumNodes)
	ConVec = read.table(confile)[,1:2]
	for ( i in 1:nrow(ConVec)){
		if(ConVec[i,1] < ConVec[i,2]+1){
			ConMatrix[ConVec[i,1],ConVec[i,2]] = 1
		}
	}
	dbs = dbscan.greater(ConMatrix,eps=eps,MinPts=MinPts,method="dist")
	
	sameclustervec = rep(0,length(GrpList))
	for( i in 1:length(GrpList)){
		ConGrpGrpVec = as.vector(unlist(GrpList[[i]]))
		NumCorePts = length(ConGrpGrpVec)
		if(! is.null(dbs$isseed)){
			if (length(which(dbs$isseed[ConGrpGrpVec])==TRUE ) == NumCorePts ){ # core points
				if (length(unique(dbs$cluster[ConGrpGrpVec])) == 1){ # same cluster
					sameclustervec[i] = 1
				}
			}
		}
	}
	
	return (list(logweight=LogWeight,sameclustervec=sameclustervec))
}
#------------------------------------------------------------------
MINPTS = 4
EPS = 1

GrpList = list()
# K cell
GrpList = c(GrpList,list(c(12,20)))
GrpList = c(GrpList,list(c(13,20)))
GrpList = c(GrpList,list(c(14,20)))
GrpList = c(GrpList,list(c(12,13,14,20)))
GrpList = c(GrpList,list(c(7,12,13,14,20)))
# GM cell
GrpList = c(GrpList,list(c(3,13)))
GrpList = c(GrpList,list(c(4,13)))
GrpList = c(GrpList,list(c(3,4,13)))
GrpList = c(GrpList,list(c(3,4,13,14,28)))

# tail
GrpList = c(GrpList,list(c(24,53)))
# NumFiles = 100

LogWeightVec = rep(0,NumFiles)
SameClusterMat = matrix(0,nrow=NumFiles,ncol=length(GrpList))
ptm <- proc.time()
for (i in 1:NumFiles){
	PtsFN = PTSFiles[i]
	ConFN = OverlapFiles[i]
	LogWeightCluster = GetWeightDensityCluster_GrpGrp(ptsfile=PtsFN,confile=ConFN,eps=EPS,MinPts=MINPTS,GrpList=GrpList)
	LogWeightVec[i] = LogWeightCluster$logweight
	SameClusterMat[i,] = LogWeightCluster$sameclustervec
}
proc.time() - ptm
WeightVec = exp(LogWeightVec - MaxLogWeight)
#------------------------------------------------------------------
OutMat = cbind(BaseNames,sprintf("%.3e",WeightVec),SameClusterMat)
FN = paste(OutDir,"/",basename(ListFN),".txt",sep="")
Header = paste("#","file","weight",paste(unlist(lapply(GrpList,function(x){paste(x,collapse="-")})), collapse="\t"),sep="\t")
write.table(file=FN,Header,quote=F,row.names=F,col.names=F,append=F)
write.table(file=FN,OutMat,quote=F,row.names=F,col.names=F,append=T,sep="\t")

# SmallestProb = min(WeightVec)/sum(WeightVec)
# Prob = rep(0,length(GrpList))
# for ( i in 1:length(GrpList)){
# 	ClusterVec = unlist(lapply(SameClusterList,function(x){x[i]}))
# 	if (sum(ClusterVec) == 0){
# 		Prob[i] = SmallestProb
# 	}else{
# 		Prob[i] = sum(WeightVec * ClusterVec) / sum(WeightVec)
# 	}
# }
#------------------------------------------------------------------
#------------------------------------------------------------------
#------------------------------------------------------------------
#------------------------------------------------------------------
