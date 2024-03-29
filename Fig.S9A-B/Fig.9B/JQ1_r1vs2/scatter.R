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


#r1_2 JQ1

JQ1_r1 <- read.table("HCT_220307_halflife_JQ1.txt", header=T)
JQ1_r2 <- read.table("HCT_240209_halflife_JQ1.txt", header=T)

###

JQ1_r1_2 <- merge(JQ1_r1, JQ1_r2, by="geneSymbol") 

heatscatter(log10(JQ1_r1_2$half_life_hr.x), log10(JQ1_r1_2$half_life_hr.y),  xlab = "log10 HL JQ1 rep1 [hr]", ylab = "log10 HL JQ1 rep2 [hr]",  xlim = c(-0.5,2), ylim = c(-0.5,2),  main="Sperman R = 0.7747716")

correlation <- cor(log10(JQ1_r1_2$half_life_hr.x), log10(JQ1_r1_2$half_life_hr.y), method = "spearman")
print(correlation)    


