#!/usr/bin/R
rm(list=ls())
graphics.off()
# ------------------------------------------------
library(scales)
library(reshape)
# ------------------------------------------------
GMCIFile="../script/tmp/GM12878.p1500_d3300.c840.con.pij.all.contactindex.txt"
KCIFile = "../script/tmp/K562.p1500_d3300.c840.con.pij.all.contactindex.txt"

GMData = read.table(GMCIFile,)
colnames(GMData) = c("NodeIndex","CI_all","CI_local","CI_long")
KData = read.table(KCIFile)
colnames(KData) = c("NodeIndex","CI_all","CI_local","CI_long")

NumNodes = nrow(GMData)
library(ggplot2)
NodeIndex = GMData[,"NodeIndex"]

# ------------------------------------------------

df_long = melt(GMData,id="NodeIndex")
p <- ggplot(data=df_long,aes(x=NodeIndex,y=value,color=variable)) +
	theme_bw() +
	geom_line() + geom_point() +
	theme(legend.position = "none") +
	labs(y = "GM Contact Index") +
	scale_colour_brewer(palette="Set1") +
	scale_x_continuous(breaks=seq(5,NumNodes,by=5)) +
	coord_cartesian(ylim=c(0,14))
dev.new(width=3, height=5)
p + facet_grid(variable ~.)
# dev.new(width=10, height=2)
# p + facet_grid(. ~ variable)

# ------------------------------------------------
df_long = melt(KData,id="NodeIndex")
p <- ggplot(data=df_long,aes(x=NodeIndex,y=value,color=variable)) +
	theme_bw() +
	geom_line() + geom_point() +
	theme(legend.position = "none") +
	labs(y = "K Contact Index") +
	scale_colour_brewer(palette="Set1") +
	scale_x_continuous(breaks=seq(5,NumNodes,by=5)) +
	coord_cartesian(ylim=c(0,14))
dev.new(width=3, height=5)
p + facet_grid(variable ~.)
# dev.new(width=10, height=2)
# p + facet_grid(. ~ variable )

# ------------------------------------------------
df = cbind(NodeIndex, (KData[,2:4] - GMData[,2:4]))
colnames(df) = c("NodeIndex","CI_all_diff","CI_local_diff","CI_long_diff")
df_long = melt(df,id="NodeIndex")
p <- ggplot(data=df_long,aes(x=NodeIndex,y=value,color=variable)) +
	theme_bw() +
	geom_line() + geom_point() +
	theme(legend.position = "none") +
	labs(y = "Contact Index Difference between K and GM") +
	scale_colour_brewer(palette="Set1") +
	scale_x_continuous(breaks=seq(5,NumNodes,by=5)) +

dev.new(width=3, height=5)
p + facet_grid(variable ~.) 