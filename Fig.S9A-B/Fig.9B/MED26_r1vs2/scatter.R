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

#r1_2 AUXIN

AUXIN_r1 <- read.table("HCT_240109_halflife_IAA.txt", header=T)
AUXIN_r2 <- read.table("HCT_240122_halflife_IAA.txt", header=T)

###

AUXIN_r1_2 <- merge(AUXIN_r1, AUXIN_r2, by="geneSymbol") 

heatscatter(log10(AUXIN_r1_2$half_life_hr.x), log10(AUXIN_r1_2$half_life_hr.y),  xlab = "log10 HL -MED26 rep1 [hr]", ylab = "log10 HL -MED26 rep2 [hr]",  xlim = c(0,1.5), ylim = c(0,1.5),  main="Sperman R = 0.4309761")

correlation <- cor(log10(AUXIN_r1_2$half_life_hr.x), log10(AUXIN_r1_2$half_life_hr.y), method = "spearman")
print(correlation)    

