
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


# Rad21 in Med26Degron (figureS10A)
dataP <- read.table("Rad21_Med26Degron.txt", sep="\t",  header=T, stringsAsFactors=F)
dataP$FC <- as.numeric(dataP$FC)

my_comparisons <- list(c("promoter", "SE"))

ggplot(dataP,aes(cat,log2(FC)))+
  geom_violin()+
  geom_boxplot(width=.1, position=position_dodge(1))+
  geom_hline(yintercept=0, lty=3, color="black")+
  geom_hline(yintercept=median(log2(dataP[dataP$cat=="promoter",]$FC)), lty=3, color="blue")+
  geom_hline(yintercept=median(log2(dataP[dataP$cat=="SE",]$FC)), lty=3, color="red")+
  labs(x='log2(FC (IAA/DMSO))', title="Rad21 signal intensity FC ")+ theme_classic()+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox", aes(label=..p.adj..),  paired = FALSE)


# Brd4 in Med26Degron 

dataP2 <- read.table("Brd4_Med26Degron.txt", sep="\t", header=T, stringsAsFactors=F)
dataP2$FC <- as.numeric(dataP2$FC)

#(figuresS10B)
ggplot(dataP2,aes(cat,log2(FC)))+
  geom_violin()+
  geom_boxplot(width=.1, position=position_dodge(1))+
  geom_hline(yintercept=0, lty=3, color="black")+
  geom_hline(yintercept=median(log2(dataP2[dataP2$cat=="promoter",]$FC)), lty=3, color="blue")+
  geom_hline(yintercept=median(log2(dataP2[dataP2$cat=="SE",]$FC)), lty=3, color="red")+
  labs(x='log2(FC (IAA/DMSO))', title="Brd4 signal intensity FC ")+ theme_classic()+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox", aes(label=..p.adj..),  paired = FALSE)

#(figureS10C)

my_comparisons <- list(c("down", "up"))


dataP3 <-subset(dataP2, cat == "promoter")
dataP3 <-subset(dataP3, catG == "up" | catG == "down")

ggplot(dataP3,aes(catG,log2(FC)))+
  geom_boxplot(outlier.shape = NA)+
  ylim(c(-2,2.5))+
  labs(x='log2(FC (IAA/DMSO))', title="Brd4 signal intensity FC at promoter") + theme_classic()+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox", aes(label=..p.adj..),  paired = FALSE)

ggplot(dataP3,aes(catG,log2(FC)))+
  geom_boxplot(outlier.shape = NA)+
  ylim(c(-1.5,1.25))+
  labs(x='log2(FC (IAA/DMSO))', title="Brd4 signal intensity FC at promoter") + theme_classic()+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox", aes(label=..p.adj..),  paired = FALSE)

