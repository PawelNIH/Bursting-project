



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

"mosaicData"

Nature <-read.csv(file = 'NatGen.csv')


Nature <- filter(Nature, pass.size == "TRUE", pass.kon == "TRUE", pass.mean == "TRUE")


ggplot(Nature, aes(x= log(mean), y= log(kon))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "spearman", size = 5) +
  geom_smooth(method='lm') 

ggplot(Nature, aes(x= log(mean), y= log(ksyn/koff))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "spearman", size = 5) +
  geom_smooth(method='lm') 



  ##
Fibro_CC <- read.csv(file = 'SKINr1_Summary_scRNA-ss-SKIN_2.csv')





Fibro_CC_SD <-filter(Fibro_CC, (Rate01MAD/Rate01Median) < 0.75,  
                          (Rate10MAD/Rate10Median) < 0.75,  
                          (EjectMAD/EjectMedian) < 0.75,    
                          Expression > 0.1)
                  



Fibro_CC_SD <- select(Fibro_CC_SD, Gene, Rate01, Eject, Rate10, Rate01Mean,Rate01Median, Rate10Median, EjectMedian, EjectMean, Rate10Mean, Decay, Expression)



Nat_CC <- merge(Nature, Fibro_CC_SD, by="Gene") 

Nat_CC_NOFILTER <- merge(Nature, Fibro_CC, by="Gene") 

ggplot(Nat_CC, aes(x= log(mean), y= log(Expression))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "pearson", size = 5) +
  geom_abline(slope=1, intercept = 0)+
  geom_smooth(method='lm')  

ggplot(Nat_CC, aes(x= log(kon), y= log(Rate01Median))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "pearson", size = 5) +
  geom_abline(slope=1, intercept = 0)+
  geom_smooth(method='lm')  


ggplot(Nat_CC, aes(x= log(kon), y= log(Rate01Median/Decay))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "pearson", size = 5) +
  geom_abline(slope=1, intercept = 0)+
  geom_smooth(method='lm')  

ggplot(Nat_CC, aes(x= log(koff), y= log(Rate10Median))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "pearson", size = 5) +
  geom_abline(slope=1, intercept = 0)+
  geom_smooth(method='lm')

ggplot(Nat_CC, aes(x= log(koff), y= log(Rate10Median/Decay))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "pearson", size = 5) +
  geom_abline(slope=1, intercept = 0)+
  geom_smooth(method='lm')


ggplot(Nat_CC, aes(x= log(ksyn), y= log(EjectMedian))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "pearson", size = 5) +
  geom_abline(slope=1, intercept = 0)+
  geom_smooth(method='lm')

ggplot(Nat_CC, aes(x= log(ksyn), y= log(EjectMedian/Decay))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "pearson", size = 5) +
  geom_abl0000ine(slope=1, intercept = 0)+
  geom_smooth(method='lm')


ggplot(Nat_CC, aes(x= log(ksyn/koff), y= log(EjectMedian/Rate10Median))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "pearson", size = 5) +
  geom_abline(slope=1, intercept = 0)+
  geom_smooth(method='lm')
##
ggplot(Nat_CC, aes(x= log(1/kon), y= log(1/(Rate01Median/Decay)))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  stat_cor(method = "pearson", size = 5) +
  geom_abline(slope=1, intercept = 0)+
  geom_smooth(method='lm')+
  geom_hex(bins = 60) +
  scale_fill_continuous(type = "viridis")


ggplot(Nat_CC, aes(x= log(ksyn/koff), y= log(EjectMedian/Rate10Median))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
 stat_cor(method = "pearson", size = 5) +
  geom_abline(slope=1, intercept = 0)+
  geom_smooth(method='lm')+
  geom_hex(bins = 60) +
  scale_fill_continuous(type = "viridis")
