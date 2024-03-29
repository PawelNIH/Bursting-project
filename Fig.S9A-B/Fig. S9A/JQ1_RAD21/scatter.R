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



DMSO_half <- read.table("HCT_merged_halflife_DMSO.txt", header=T)
AUXIN_half <- read.table("HCT_merged_halflife_IAA.txt", header=T)
JQ1_half <- read.table("HCT_merged_halflife_JQ1.txt", header=T)



###

DMSO_RAD21 <- merge(DMSO_half, AUXIN_half, by="geneSymbol") 


log10_half_life_DMSO <- log10(DMSO_RAD21$half_life_hr.x)
log10_half_life_AUXIN <- log10(DMSO_RAD21$half_life_hr.y)

heatscatter(log10_half_life_DMSO, log10_half_life_AUXIN, xlab = "log10 HL DMSO [hr]", ylab = "log10 HL -RAD21 [hr]",  xlim = c(-1,2.5), ylim = c(-1,2.5),  main="Sperman R = 0.8273703")

correlation <- cor(log10_half_life_DMSO, log10_half_life_AUXIN, method = "spearman")
print(correlation)    

###

DMSO_JQ1 <- merge(DMSO_half, JQ1_half, by="geneSymbol") 


log10_half_life_DMSO <- log10(DMSO_JQ1$half_life_hr.x)
log10_half_life_JQ1 <- log10(DMSO_JQ1$half_life_hr.y)

heatscatter(log10_half_life_DMSO, log10_half_life_JQ1, xlab = "log10 HL DMSO [hr]", ylab = "log10 HL JQ1 [hr]",  xlim = c(-1,2.5), ylim = c(-1,2.5),  main="Sperman R = 0.7846899")

correlation <- cor(log10_half_life_DMSO, log10_half_life_JQ1, method = "spearman")
print(correlation)    


