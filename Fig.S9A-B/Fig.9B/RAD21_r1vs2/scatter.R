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

AUXIN_r1 <- read.table("HCT_220307_halflife_IAA.txt", header=T)
AUXIN_r2 <- read.table("HCT_240209_halflife_IAA.txt", header=T)

###

AUXIN_r1_2 <- merge(AUXIN_r1, AUXIN_r2, by="geneSymbol") 

heatscatter(log10(AUXIN_r1_2$half_life_hr.x), log10(AUXIN_r1_2$half_life_hr.y),  xlab = "log10 HL -RAD21 rep1 [hr]", ylab = "log10 HL -RAD21 rep2 [hr]",  xlim = c(-0.5,2), ylim = c(-0.5,2),  main="Sperman R = 0.7566393")

correlation <- cor(log10(AUXIN_r1_2$half_life_hr.x), log10(AUXIN_r1_2$half_life_hr.y), method = "spearman")
print(correlation)    



