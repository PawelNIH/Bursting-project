
library('readxl')
library('ggplot2')
library('ggrepel')
library('plyr')
library('stringr')
library('knitr')
library('data.table')
library('ggpubr')
library("RColorBrewer")
library('dplyr')
library("gghighlight")
library("LSD")
library("ggstance")
library("mosaicData")

#JQ1
PI2 <- read.table("Pausing_index_JQ1.txt", sep="\t", header=T, stringsAsFactors=F)
PI2$type3 <- PI2$type2
PI2$type3 <- ifelse( PI2$type2 == "down", "down", "other")

my_comparisons <- list(c("down", "other"))


ggplot(PI2, aes(type3,LFC, color=type3))+ geom_violin()+
	geom_boxplot(width=.1)+
	ylim(c(-6,7))+
	geom_hline(yintercept=median(PI2[PI2$type2=="down",]$LFC, na.rm=T), lty=3, color="black")+
	labs(title="Pausing index log2FC") + theme_classic()+   
  scale_color_brewer(palette="Set1") +
  stat_compare_means(comparisons=my_comparisons, method = "wilcox", aes(label=..p.adj..),  paired = FALSE)

table(PI2$type3)

#Med26Degron
PI1 <- read.table("Pausing_index_Med26Degron.txt", sep="\t", header=T, stringsAsFactors=F)

PI1$type3 <- PI1$type2
PI1$type3 <- ifelse( PI1$type2 == "down", "down", "other")

ggplot(PI1, aes(type3,LFC, color=type3))+geom_violin()+
	geom_boxplot(width=.1)+
	ylim(c(-6,7))+
	geom_hline(yintercept=median(PI1[PI1$type2=="down",]$LFC, na.rm=T), lty=3, color="black")+
	labs(title="Pausing index log2FC") + theme_classic()+   
  scale_color_brewer(palette="Set1")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox", aes(label=..p.adj..),   paired = FALSE,
)

table(PI1$type3)

