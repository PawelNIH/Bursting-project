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


##Load scRNA-Seq data

all_DMSO <- read.csv(file = 'all_DMSO_Summary_scRNA-ss-HCT116_2.csv')

all_burst_DMSO <- read.csv(file = 'all_DMSO_burst_scRNA-ss-HCT116_HCTDMSO2h_2.csv')


all_DMSO <- merge(all_DMSO, all_burst_DMSO, by="Gene") 

##


all_AUXIN <- read.csv(file = 'all_AUXIN_Summary_scRNA-ss-HCT116_2.csv')

all_burst_AUXIN <- read.csv(file = 'all_AUXIN_burst_scRNA-ss-HCT116_HCTAUXIN2h_2.csv')


all_AUXIN <- merge(all_AUXIN, all_burst_AUXIN, by="Gene") 

##
sc <-  dplyr::bind_rows(all_DMSO, all_AUXIN)

##QC filtering
sc <-filter(sc, (Rate01MAD/Rate01Median) < 0.75,  
            (Rate10MAD/Rate10Median) < 0.75,  
            (EjectMAD/EjectMedian) < 0.75,  
            (BurstMAD/BurstMedian)< 0.75,  
            Expression > 0.01)

sc <- transform(sc, log10_Rate01Median = log10(Rate01Median))

sc <- transform(sc, log10_Rate01MADplus = log10(Rate01Median + Rate01MAD))

sc <- transform(sc, log10_Rate01MADminus = log10(Rate01Median - Rate01MAD))



#

sc <- transform(sc, log10_Rate10Median = log10(Rate10Median))

sc <- transform(sc, log10_Rate10MADplus = log10(Rate10Median + Rate10MAD))

sc <- transform(sc, log10_Rate10MADminus = log10(Rate10Median - Rate10MAD))


#

sc <- transform(sc, log10_EjectMedian = log10(EjectMedian))

sc <- transform(sc, log10_EjectMADplus = log10(EjectMedian + EjectMAD))

sc <- transform(sc, log10_EjectMADminus = log10(EjectMedian - EjectMAD))

#
#

sc <- transform(sc, log10_BurstMedian = log10(BurstMedian))

sc <- transform(sc, log10_BurstMADplus = log10(BurstMedian + BurstMAD))

sc <- transform(sc, log10_BurstMADminus = log10(BurstMedian - BurstMAD))


##
sc <- select(sc, Gene, log10_Rate01Median, log10_Rate01MADplus, log10_Rate01MADminus, log10_Rate10Median, log10_Rate10MADplus, log10_Rate10MADminus,  log10_EjectMedian, log10_EjectMADplus, log10_EjectMADminus, log10_BurstMedian, log10_BurstMADplus, log10_BurstMADminus, Condition, Expression)


sc <- filter(sc, Gene %in% c("MYC", "SOX9", "ERRFI1", "HMMR", "KIF2C", "PPM1G", "PSMD14"))


## load smFISH data
FISH <- read.csv(file = 'Summary_FISH-ss-HCT116_2.csv')

FISH_DMSO <- filter(FISH, Condition == "DMSO")

FISH_burst_DMSO <- read.csv(file = 'burst_FISH-ss-HCT116_DMSO_2.csv')

FISH_DMSO <- merge(FISH_DMSO, FISH_burst_DMSO, by="Gene") 

##

FISH_AUXIN <- filter(FISH, Condition == "AUXIN")


FISH_burst_AUXIN <- read.csv(file = 'burst_FISH-ss-HCT116_IAA_2.csv')

FISH_AUXIN <- merge(FISH_AUXIN, FISH_burst_AUXIN, by="Gene") 


FISH <-  dplyr::bind_rows(FISH_DMSO, FISH_AUXIN)

FISH <- transform(FISH, log10_Rate01Median = log10(Rate01Median))

FISH <- transform(FISH, log10_Rate01MADplus = log10(Rate01Median + Rate01MAD))

FISH <- transform(FISH, log10_Rate01MADminus = log10(Rate01Median - Rate01MAD))

##

FISH <- transform(FISH, log10_Rate10Median = log10(Rate10Median))

FISH <- transform(FISH, log10_Rate10MADplus = log10(Rate10Median + Rate10MAD))

FISH <- transform(FISH, log10_Rate10MADminus = log10(Rate10Median - Rate10MAD))

##

#

FISH <- transform(FISH, log10_EjectMedian = log10(EjectMedian))

FISH <- transform(FISH, log10_EjectMADplus = log10(EjectMedian + EjectMAD))

FISH <- transform(FISH, log10_EjectMADminus = log10(EjectMedian - EjectMAD))

##

FISH <- transform(FISH, log10_BurstMedian = log10(BurstMedian))

FISH <- transform(FISH, log10_BurstMADplus = log10(BurstMedian + BurstMAD))

FISH <- transform(FISH, log10_BurstMADminus = log10(BurstMedian - BurstMAD))


##
FISH <- select(FISH, Gene, log10_Rate01Median, log10_Rate01MADplus, log10_Rate01MADminus, log10_Rate10Median, log10_Rate10MADplus, log10_Rate10MADminus,  log10_EjectMedian, log10_EjectMADplus, log10_EjectMADminus, log10_BurstMedian, log10_BurstMADplus, log10_BurstMADminus, Condition, Expression)

FISH <- filter(FISH, Gene %in% c("MYC", "SOX9", "ERRFI1", "HMMR", "KIF2C", "PPM1G", "PSMD14"))


sc_FISH <-  bind_cols(sc, FISH)

##

library("lmodel2")
library("ggpmisc")


lmodel2(formula = log10_Rate01Median...2 ~ log10_Rate01Median...17, data = sc_FISH,
        range.y = "interval", range.x = "interval", nperm = 99)

lmodel2(formula = log10_Rate10Median...5 ~ log10_Rate10Median...20, data = sc_FISH,
        range.y = "interval", range.x = "interval", nperm = 99)

lmodel2(formula = log10_EjectMedian...8 ~ log10_EjectMedian...23, data = sc_FISH,
        range.y = "interval", range.x = "interval", nperm = 99)

lmodel2(formula = log10_BurstMedian...11 ~ log10_BurstMedian...26, data = sc_FISH,
        range.y = "interval", range.x = "interval", nperm = 99)


####################################### color by  scRNA expression


ggplot(sc_FISH, aes(y= log10_Rate01Median...2, x= log10_Rate01Median...17, color = Expression...15)) +  geom_point(size = 2) +
  geom_smooth(method='lm', se = F) +theme_classic()+ 
  geom_errorbar(aes(ymin = log10_Rate01MADminus...4,ymax = log10_Rate01MADplus...3)) + 
  geom_errorbarh(aes(xmin = log10_Rate01MADminus...19,xmax = log10_Rate01MADplus...18))+scale_color_gradientn(colours = rainbow(5))+
  stat_cor(method = "pearson", aes(label = ..r.label..))


#



ggplot(sc_FISH, aes(y= log10_Rate10Median...5, x= log10_Rate10Median...20, color = Expression...15)) +  geom_point(size = 2) +
  geom_smooth(method='lm', se = F) +theme_classic()+ 
  geom_errorbar(aes(ymin = log10_Rate10MADminus...7,ymax = log10_Rate10MADplus...6)) + 
  geom_errorbarh(aes(xmin = log10_Rate10MADminus...22,xmax = log10_Rate10MADplus...21))+scale_color_gradientn(colours = rainbow(5))+
  stat_cor(method = "pearson", aes(label = ..r.label..))



#


ggplot(sc_FISH, aes(y= log10_EjectMedian...8, x= log10_EjectMedian...23, color = Expression...15)) +  geom_point(size = 2) +
  geom_smooth(method='lm', se = F) +theme_classic()+ 
  geom_errorbar(aes(ymin = log10_EjectMADminus...10,ymax = log10_EjectMADplus...9)) + 
  geom_errorbarh(aes(xmin = log10_EjectMADminus...25,xmax = log10_EjectMADplus...24))+scale_color_gradientn(colours = rainbow(5))+
  stat_cor(method = "pearson", aes(label = ..r.label..))





#

ggplot(sc_FISH, aes(y= log10_BurstMedian...11, x= log10_BurstMedian...26, color = Expression...15)) +  geom_point(size = 2) +
  geom_smooth(method='lm', se = F) +theme_classic()+ 
  geom_errorbar(aes(ymin = log10_BurstMADminus...13,ymax = log10_BurstMADplus...12)) + 
  geom_errorbarh(aes(xmin = log10_BurstMADminus...28,xmax = log10_BurstMADplus...27))+scale_color_gradientn(colours = rainbow(5))+
  stat_cor(method = "pearson", aes(label = ..r.label..))


