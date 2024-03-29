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

##


all_DMSO <- read.csv(file = "DMSOr1_Summary_scRNA-ss-HCT116_2.csv")

all_burst_DMSO <- read.csv(file = 'DMSOr1_burst_scRNA-ss-HCT116_MED26DDMSOr1_2.csv')


all_DMSO <- merge(all_DMSO, all_burst_DMSO, by="Gene") 

##


all_DMSOr2 <- read.csv(file = "DMSOr2_Dr2_Summary_scRNA-ss-HCT116_2.csv")

all_burst_DMSOr2 <- read.csv(file = "DMSOr2_burst_scRNA-ss-HCT116_MED26DDMSOr2_2.csv")


all_DMSOr2 <- merge(all_DMSOr2, all_burst_DMSOr2, by="Gene") 

##

allr12 <- merge(all_DMSO, all_DMSOr2, by="Gene") 
##

allr12_filter <-filter(allr12, (Rate01MAD.x/Rate01Median.x) < 0.75,  (Rate01MAD.y/Rate01Median.y) < 0.75,  
                            (Rate10MAD.x/Rate10Median.x) < 0.75,  (Rate10MAD.y/Rate10Median.y) < 0.75, 
                            (EjectMAD.x/EjectMedian.x) < 0.75,  (EjectMAD.y/EjectMedian.y) < 0.75,   
                            (BurstMAD.x/BurstMedian.x) < 0.75, (BurstMAD.y/BurstMedian.y) < 0.75,  
                            Expression.x > 0.01, Expression.y > 0.01)

ggplot(allr12_filter, aes(x= log10(Rate01Median.x), y= log10(Rate01Median.y))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "pearson", size = 5) +
  geom_abline(slope=1, intercept = 0)+
  geom_smooth(method='lm')  


ggplot(allr12_filter, aes(x= log10(Rate10Median.x), y= log10(Rate10Median.y))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "pearson", size = 5) +
  geom_abline(slope=1, intercept = 0)+
  geom_smooth(method='lm')  


ggplot(allr12_filter, aes(x= log10(EjectMedian.x), y= log10(EjectMedian.y))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "pearson", size = 5) +
  geom_abline(slope=1, intercept = 0)+
  geom_smooth(method='lm')  



ggplot(allr12_filter, aes(x= log10(BurstMedian.x), y= log10(BurstMedian.y))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "pearson", size = 5) +
  geom_abline(slope=1, intercept = 0)+
  geom_smooth(method='lm')  

ggplot(allr12_filter, aes(x= log10(Expression.x), y= log10(Expression.y))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "pearson", size = 5) +
  geom_abline(slope=1, intercept = 0)+
  geom_smooth(method='lm') 



##low
EXP_allr12_filter <- filter(allr12_filter,  Expression.x < 1,  Expression.y < 1)

##middle

EXP_allr12_filter <- filter(allr12_filter, Expression.x > 1, Expression.x < 5, Expression.y > 1,  Expression.y < 5)

##high

EXP_allr12_filter <- filter(allr12_filter, Expression.x > 5, Expression.y > 5)

##
ggplot(EXP_allr12_filter, aes(x= log10(Rate01Median.x), y= log10(Rate01Median.y))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "pearson", size = 5) +
  geom_abline(slope=1, intercept = 0)+
  geom_smooth(method='lm')  


ggplot(EXP_allr12_filter, aes(x= log10(Rate10Median.x), y= log10(Rate10Median.y))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "pearson", size = 5) +
  geom_abline(slope=1, intercept = 0)+
  geom_smooth(method='lm')  


ggplot(EXP_allr12_filter, aes(x= log10(EjectMedian.x), y= log10(EjectMedian.y))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "pearson", size = 5) +
  geom_abline(slope=1, intercept = 0)+
  geom_smooth(method='lm')  


ggplot(EXP_allr12_filter, aes(x= log10(BurstMedian.x), y= log10(BurstMedian.y))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "pearson", size = 5) +
  geom_abline(slope=1, intercept = 0)+
  geom_smooth(method='lm')  

ggplot(EXP_allr12_filter, aes(x= log10(Expression.x), y= log10(Expression.y))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "pearson", size = 5) +
  geom_abline(slope=1, intercept = 0)+
  geom_smooth(method='lm')  

