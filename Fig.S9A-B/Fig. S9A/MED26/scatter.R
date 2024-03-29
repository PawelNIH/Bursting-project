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

#



DMSO_half <- read.table("HCT_Med26Degron_halflife_DMSO.txt", header=T)
AUXIN_half <- read.table("HCT_Med26Degron_halflife_IAA.txt", header=T)



###

DMSO_MED26 <- merge(DMSO_half, AUXIN_half, by="geneSymbol") 



log10_half_life_DMSO <- log10(DMSO_MED26$half_life_hr.x)
log10_half_life_AUXIN <- log10(DMSO_MED26$half_life_hr.y)

heatscatter(log10_half_life_DMSO, log10_half_life_AUXIN, xlab = "log10 HL DMSO [hr]", ylab = "log10 HL -MED26 [hr]",  xlim = c(-0.8,2), ylim = c(-0.8,2),  main="Sperman R = 0.5074518")


correlation <- cor(log10_half_life_DMSO, log10_half_life_AUXIN, method = "spearman")
print(correlation)    

