#!/usr/bin/R

rm(list=ls())
PtsFile = "GM.pts"
PtsFile = "0089.pts"
PtsFile = "/tmp/AA/0001.pts"
PtsFile = "/tmp/AA_0001.pts"
ConIndDstFile = "K562.p1500_d3300.q5.a5.sis.node.dst"
# ConIndDstFile = "GM12878.p1500_d3300.q5.a5.sis.node.dst"

Pts = read.table(PtsFile)
ConIndDst = read.table(ConIndDstFile)
Con.Vec = ConIndDst[,1:2]

Distance <- function(p1,p2){
	Diff = p1[1:3] - p2[1:3]
	return(sqrt(sum(Diff*Diff)))
}

Dist.Vec = apply(Con.Vec,1,function(x){Distance(Pts[x[1],],Pts[x[2],])})

Diff <- Dist.Vec - ConIndDst[,3]
mean(Diff)