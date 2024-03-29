
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



my_comparisons <- list(c("down", "others"))

##Rad21 degron

data <- read.table("TBP_Rad21Degron.txt", header = TRUE, stringsAsFactors = FALSE)
data$FC <- as.numeric(data$FC)

ggplot(data,aes(catO,log2(FC)))+
  geom_violin()+
  geom_boxplot(width=.1, position=position_dodge(1))+
  geom_hline(yintercept=0, lty=3, color="black")+
  geom_hline(yintercept=median(log2(data[data$catO=="down",]$FC)), lty=3, color="blue")+
  ylim(-4,3)+
  labs(x='log2(FC (Auxin/DMSO))', title="TBP signal intensity FC ")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox", aes(label=..p.adj..),  paired = FALSE)


## JQ1 treatment
data <- read.table("TBP_JQ1.txt", header = TRUE, stringsAsFactors = FALSE)
data$FC <- as.numeric(data$FC)

ggplot(data,aes(catO,log2(FC)))+
  geom_violin()+
  geom_boxplot(width=.1, position=position_dodge(1))+
  geom_hline(yintercept=0, lty=3, color="black")+
  geom_hline(yintercept=median(log2(data[data$catO=="down",]$FC)), lty=3, color="blue")+
  ylim(-4,3)+
  labs(x='log2(FC (JQ1/DMSO))', title="TBP signal intensity FC ")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox", aes(label=..p.adj..),  paired = FALSE)




##Med26 degron
data <- read.table("TBP_Med26Degron.txt", header = TRUE, stringsAsFactors = FALSE)
data$FC <- as.numeric(data$FC)

ggplot(data,aes(catO,log2(FC)))+
  geom_violin()+
  geom_boxplot(width=.1, position=position_dodge(1))+
  geom_hline(yintercept=0, lty=3, color="black")+
  geom_hline(yintercept=median(log2(data[data$catO=="down",]$FC)), lty=3, color="blue")+
  ylim(-4,4)+
  labs(x='log2(FC (IAA/DMSO))', title="TBP signal intensity FC ")+
  stat_compare_means(comparisons=my_comparisons, method = "wilcox", aes(label=..p.adj..),  paired = FALSE)



