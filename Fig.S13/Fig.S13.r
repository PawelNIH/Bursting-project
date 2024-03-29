library("ggplot2")
library("plyr")
library("grid")
library('Hmisc')
library("gplots")
library("dynamicTreeCut")
library("RColorBrewer")
library("boot")
library('DESeq2')
library("lattice")
library("latticeExtra")
library('hexbin')
library('reshape')
library('LSD')
library('psych')

annot <- data.frame(lab= c("-5kb", "Start", "End", "+5kb"), mpos = c(0,49,149,199))
lines <- data.frame( x = c(49, 149))

## JQ1
rpkmA <- read.table("PolII_HCTJQ1_data.txt", sep="\t", header=T, stringsAsFactors=F)

#other category
temp <- subset(rpkmA, gt=="OJQ1" | gt=="ODMSO")
ggplot(aes(pos,Mn, colour=gt), data=temp) +
  geom_line()+
  geom_vline(aes(xintercept = x), data = lines, linetype = 3) +
  scale_x_continuous(breaks = annot$mpos, labels = annot$lab)+
  labs(x = "Position", y= "mean rpkm", title = "polII composite on each category for Other genes")

#down category
temp <- subset(rpkmA, gt=="DJQ1" | gt=="DDMSO")
ggplot(aes(pos,Mn, colour=gt), data=temp) +
  geom_line()+
  geom_vline(aes(xintercept = x), data = lines, linetype = 3) +
  scale_x_continuous(breaks = annot$mpos, labels = annot$lab)+
  ylim(c(0,20))+
  theme_classic()+
  labs(x = "Position", y= "mean rpkm", title = "polII composite on each category for Down Sig genes")

##Rad21 Degron
rpkmA <- read.table("PolII_HCTRad21Degron_data.txt", sep="\t", header=T, stringsAsFactors=F)

#other category
temp <- subset(rpkmA, gt=="OIAA" | gt=="ODMSO")
ggplot(aes(pos,Mn, colour=gt), data=temp) +
  geom_line()+
  geom_vline(aes(xintercept = x), data = lines, linetype = 3) +
  scale_x_continuous(breaks = annot$mpos, labels = annot$lab)+
  labs(x = "Position", y= "mean rpkm", title = "polII composite on each category for Other genes")

#down category
temp <- subset(rpkmA, gt=="DIAA" | gt=="DDMSO")
ggplot(aes(pos,Mn, colour=gt), data=temp) +
  geom_line()+
  geom_vline(aes(xintercept = x), data = lines, linetype = 3) +
  scale_x_continuous(breaks = annot$mpos, labels = annot$lab)+
  ylim(c(0,20))+
  theme_classic()+
  labs(x = "Position", y= "mean rpkm", title = "polII composite on each category for Down Sig genes")

## Med26 Degron
rpkmA <- read.table("PolII_HCTMed26Degron_data.txt", sep="\t", header=T, stringsAsFactors=F)

#other category
temp <- subset(rpkmA, gt=="OIAA" | gt=="ODMSO")
ggplot(aes(pos,Mn, colour=gt), data=temp) +
  geom_line()+
  geom_vline(aes(xintercept = x), data = lines, linetype = 3) +
  scale_x_continuous(breaks = annot$mpos, labels = annot$lab)+
  labs(x = "Position", y= "mean rpkm", title = "polII composite on each category for Other genes")

#down category
temp <- subset(rpkmA, gt=="DIAA" | gt=="DDMSO")
ggplot(aes(pos,Mn, colour=gt), data=temp) +
  geom_line()+
  geom_vline(aes(xintercept = x), data = lines, linetype = 3) +
  scale_x_continuous(breaks = annot$mpos, labels = annot$lab)+
  ylim(c(0,20))+
  theme_classic()+
  labs(x = "Position", y= "mean rpkm", title = "polII composite on each category for Down Sig genes")
