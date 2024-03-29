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


##Low_vs_high_rep1

lowDMSOr1_summary <- read.csv(file = 'all_DMSO_r1_Summary_scRNA-ss-HCT116_2.csv')

lowDMSOr1_burst <- read.csv(file = 'all_DMSO_r1_burst_scRNA-ss-HCT116_MED26DDMSOr1_2.csv')

lowDMSOr1 <- merge(lowDMSOr1_summary, lowDMSOr1_burst, by="Gene") 


lowDMSOr1_MAD <-filter(lowDMSOr1, (Rate01MAD/Rate01Median) < 0.75, 
                                   (Rate10MAD/Rate10Median) < 0.75,  
                                   (EjectMAD/EjectMedian) < 0.75,  
                       BurstMAD/BurstMedian < 0.75, 
                      Expression > 0.01)



highDMSOr1_summary <- read.csv(file = "high_all_DMSOr1_Summary_scRNA-ss-HCT116_2.csv")

highDMSOr1_burst <- read.csv(file = " high_all_DMSOr1_burst_MED26DDMSOr1_2.csv")

highDMSOr1 <- merge(highDMSOr1_summary, highDMSOr1_burst, by="Gene") 


highDMSOr1_MAD <-filter(highDMSOr1, (Rate01MAD/Rate01Median) < 0.75, 
                       (Rate10MAD/Rate10Median) < 0.75,  
                       (EjectMAD/EjectMedian) < 0.75,  
                       BurstMAD/BurstMedian < 0.75, 
                       Expression > 0.01)


##
DMSOr1_filtered <- merge(lowDMSOr1_MAD, highDMSOr1_MAD, by="Gene") 


##

ggplot(DMSOr1_filtered, aes(x= log10(Rate01Median.x), y= log10(Rate01Median.y))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "pearson", size = 5) +
  geom_abline(slope=1, intercept = 0, color = "red", size = 2)+
  geom_bin2d(bins = 70) + xlab("ON rate lower depth (log10)")+ ylab("ON rate higher depth (log10)")


ggplot(DMSOr1_filtered, aes(x= log10(BurstMedian.x), y= log10(BurstMedian.y))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "pearson", size = 5) +
  geom_abline(slope=1, intercept = 0, color = "red", size = 2)+
  geom_bin2d(bins = 70)+ xlab("Burst Size lower depth (log10")+ ylab("Burst Size higher depth (log10)")


ggplot(DMSOr1_filtered, aes(x= log10((BurstMedian.x/42760.3674)*79327.7048), y= log10(BurstMedian.y))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "pearson", size = 5) +
  geom_abline(slope=1, intercept = 0, color = "red", size = 2)+
  geom_bin2d(bins = 70)+ xlab("Normalized Burst Size lower depth (log10)")+ ylab("Burst Size higher depth (log10)")


