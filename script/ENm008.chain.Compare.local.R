#!/usr/bin/R
# Compare.local.R
rm(list=ls())
graphics.off()
# ------------------------------------------------
library(scales)
library(reshape)
# ------------------------------------------------
Node1 = 3
Node2 = 35

GMFile = paste("../script/local/GM12878.p1500_d3300.c840.r270.ovl.local.",Node1,"-",Node2,".txt",sep="")
KFile = paste("../script/local/K562.p1500_d3300.c840.r270.ovl.local.",Node1,"-",Node2,".txt",sep="")


GMData = read.table(GMFile)
colnames(GMData) = c("NodeIndex","Triplet",paste("K-",Node1,"-",Node2,sep=""),paste(Node1,"-",Node2,"-K",sep=""),"Intercept")
KData = read.table(KFile)
colnames(KData) = c("NodeIndex","Triplet",paste("K-",Node1,"-",Node2,sep=""),paste(Node1,"-",Node2,"-K",sep=""),"Intercept")

NumNodes = nrow(GMData)
library(ggplot2)
NodeIndex = GMData[,"NodeIndex"]

# ------------------------------------------------
# ------------------------------------------------
df_long = melt(GMData,id="NodeIndex")
vline.data1 <- data.frame(variable = levels(df_long$variable), vl=c(Node1,Node1,Node2,Node1)) 
vline.data2 <- data.frame(variable = levels(df_long$variable), vl=c(Node2,Node1,Node2,Node2)) 

p1 <- ggplot(data=df_long,aes(x=NodeIndex,y=value,color=variable)) +
	theme_bw() +
	geom_line() + geom_point() +
	theme(legend.position = "none") +
	labs(title="GM12878", x= paste("Node Index"),y = paste("Node (",Node1 , "-", Node2,") Contact Propensity")) + 
	scale_colour_brewer(palette="Set1") +
	scale_x_continuous(breaks=seq(5,NumNodes,by=5)) +
	scale_y_continuous(breaks=seq(0,0.6,by=0.2)) +
	coord_cartesian(ylim=c(0,0.8))
	
# dev.new(width=3, height=5)
p1 <- p1 + facet_grid(variable ~.) +
geom_vline(aes(xintercept = vl), vline.data1,linetype="dashed",size=.5,colour="grey40") +
geom_vline(aes(xintercept = vl), vline.data2,linetype="dashed",size=.5,colour="grey40") 

# dev.new(width=10, height=2)
# p + facet_grid(. ~ variable)

# ------------------------------------------------

df_long = melt(KData,id="NodeIndex")
p2 <- ggplot(data=df_long,aes(x=NodeIndex,y=value,color=variable)) +
	theme_bw() +
	geom_line() + geom_point() +
	theme(legend.position = "none") +
	labs(title="K562", x= paste("Node Index"),y = paste("Node (",Node1 , "-", Node2,") Contact Propensity")) + 
	scale_colour_brewer(palette="Set1") +
	scale_x_continuous(breaks=seq(5,NumNodes,by=5)) +
	scale_y_continuous(breaks=seq(0,0.6,by=0.2)) +
	coord_cartesian(ylim=c(0,0.8))
	
# dev.new(width=3, height=5)
p2 <- p2 + facet_grid(variable ~.) +
geom_vline(aes(xintercept = vl), vline.data1,linetype="dashed",size=.5,colour="grey40") +
geom_vline(aes(xintercept = vl), vline.data2,linetype="dashed",size=.5,colour="grey40") 


# ------------------------------------------------
df = cbind(NodeIndex, (KData[,2:5] - GMData[,2:5]))
df_long = melt(df,id="NodeIndex")

p3 <- ggplot(data=df_long,aes(x=NodeIndex,y=value,color=variable)) +
	theme_bw() +
	geom_line() + geom_point() +
	theme(legend.position = "none") +
	labs(title="K562 - GM12878", x= paste("Node Index"),y = paste("Node (",Node1 , "-", Node2,") Contact Propensity")) + 
	scale_colour_brewer(palette="Set1") +
	scale_x_continuous(breaks=seq(5,NumNodes,by=5))  +
	scale_y_continuous(breaks=seq(-0.1,0.1,by=0.1)) +
	coord_cartesian(ylim=c(-0.2,0.2))
	
# dev.new(width=3, height=5)
p3 <- p3 + facet_grid(variable ~.) +
geom_vline(aes(xintercept = vl), vline.data1,linetype="dashed",size=.5,colour="grey40") +
geom_vline(aes(xintercept = vl), vline.data2,linetype="dashed",size=.5,colour="grey40") 

library(gridExtra)
dev.new(width=12,height=7)
grid.arrange(p1,p2,p3,ncol=3)
