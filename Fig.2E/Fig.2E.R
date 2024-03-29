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


####SMAD3
SMAD3_raw <- read.table("SMAD3_HCTDMSOr1.txt", header=F)

SMAD3_pred <- read.table("SMAD3pred.txt", header=F)

##Plot
plot(c(0:11),SMAD3_raw$V1,type="h", col = "gray", lwd = 4, xlim=c(0, 10))

lines(c(0:19),SMAD3_pred$V*1972,type="l", col = "red")

##Bursting

SMAD3 <- read.csv(file = 'SMAD3.csv')

SMAD3 <- transform(SMAD3, off_time = 1/Rate01Median)

SMAD3 <- transform(SMAD3, off_time_MAD = (1/Rate01Median^2*Rate01MAD))

SMAD3_off_time <-SMAD3$off_time

SMAD3_off_time_MAD <-SMAD3$off_time_MAD

SMAD3_burst_size <-SMAD3$BurstMedian

SMAD3_burst_MAD <-SMAD3$BurstMAD

####POLR2K

POLR2K_raw <- read.table("POLR2K_HCTDMSOr1.txt", header=F)

POLR2K_pred <- read.table("POLR2Kpred.txt", header=F)

##Plot

POLR2K_raw <- POLR2K_raw[c(0:30),]

plot(c(0:29),POLR2K_raw$V1,type="h", col = "gray", lwd = 4, xlim=c(0, 22))

lines(c(0:29),POLR2K_pred$V1*1972,type="l", col = "red")

###Bursting

POLR2K <- read.csv(file = 'POLR2K.csv')

POLR2K <- transform(POLR2K, off_time = 1/Rate01Median)

POLR2K <- transform(POLR2K, off_time_MAD = (1/Rate01Median^2*Rate01MAD))

POLR2K_off_time <-POLR2K$off_time

POLR2K_off_time_MAD <-POLR2K$off_time_MAD

POLR2K_burst_size <-POLR2K$BurstMedian

POLR2K_burst_MAD <-POLR2K$BurstMAD
