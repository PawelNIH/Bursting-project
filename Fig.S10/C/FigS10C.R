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

##Transcriptional bursting analyses 


all_DMSO <- read.csv(file = 'DMSO_Summary_scRNA-ss-ESC_2.csv')

all_burst_DMSO <- read.csv(file = 'DMSO_burst_scRNA-ss-ESC_ESCRDMSO_2.csv')


all_DMSO <- merge(all_DMSO, all_burst_DMSO, by="Gene") 


##


all_AUXIN <- read.csv(file = 'AUXIN_Summary_scRNA-ss-ESC_2.csv')

all_burst_AUXIN <- read.csv(file = 'AUXIN_burst_scRNA-ss-ESC_ESCRAUXIN_2.csv')


all_AUXIN <- merge(all_AUXIN, all_burst_AUXIN, by="Gene") 


##
DMSO_AUXIN_horizontal <- merge(all_DMSO, all_AUXIN, by="Gene") 


##QC

DMSO_AUXIN_horizontal_filter <-filter(DMSO_AUXIN_horizontal, (Rate01MAD.x/Rate01Median.x) < 0.75,  (Rate01MAD.y/Rate01Median.y) < 0.75 ,  
                                      (Rate10MAD.x/Rate10Median.x) < 0.75,  (Rate10MAD.y/Rate10Median.y) < 0.75 , 
                                      (EjectMAD.x/EjectMedian.x) < 0.75,  (EjectMAD.y/EjectMedian.y) < 0.75, 
                                      (BurstMAD.x/BurstMedian.x) < 0.75, (BurstMAD.y/BurstMedian.y) < 0.75,    Expression.x > 0.01,  Expression.y > 0.01)



##
ggplot(DMSO_AUXIN_horizontal_filter, aes(x = log10((1/Rate01Median.y)/(1/Rate01Median.x)), y = log10(Expression.y/Expression.x))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) +
  xlab("log10 (Treatment/Ctrl) OFF time") + ylab("log10 (Treatment/Ctrl) Expression") +
  stat_cor(method = "pearson", size = 5) +
  geom_smooth(method='lm')



ggplot(DMSO_AUXIN_horizontal_filter, aes(x = log10(BurstMedian.y/BurstMedian.x), y = log10(Expression.y/Expression.x), )) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + xlab("log10 (Treatment/Ctrl) Burst Size") + ylab("log10 (Treatment/Ctrl) Expression") +
  stat_cor(method = "pearson", size = 5) +
  geom_smooth(method='lm')

