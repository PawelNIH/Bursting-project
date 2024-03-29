
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


all_DMSO <- read.csv(file = 'all_Dr2_Summary_scRNA-ss-HCT116_2.csv')

all_burst_DMSO <- read.csv(file = 'all_Dr2_burst_scRNA-ss-HCT116_MED26DDMSOr2_2.csv')


all_DMSO <- merge(all_DMSO, all_burst_DMSO, by="Gene") 


##


all_AUXIN <- read.csv(file = 'all_Ir2_Summary_scRNA-ss-HCT116_2.csv')

all_burst_AUXIN <- read.csv(file = 'all_Ir2_burst_scRNA-ss-HCT116_MED26DIAAr2_2.csv')


all_AUXIN <- merge(all_AUXIN, all_burst_AUXIN, by="Gene") 


##
DMSO_AUXIN <-  dplyr::bind_rows(all_DMSO, all_AUXIN)

##
DMSO_AUXIN <- transform(DMSO_AUXIN, off_time = (1/Rate01Median))


## rnaseq IAA bulk categories T120
rnaseq <- read.table("HCT_Med26Degron_fdroutput_list", sep=",", header=T)
down_120 <- subset(rnaseq,  ud == "down")

DMSO_AUXIN$cat_IAA_120_bulk <- ifelse(DMSO_AUXIN$Gene %in% down_120$refseq_name, "down", "others")

##
DMSO_AUXIN_horizontal <- merge(all_DMSO, all_AUXIN, by="Gene") 

DMSO_AUXIN_horizontal$cat_IAA_120_bulk <- ifelse(DMSO_AUXIN_horizontal$Gene %in% down_120$refseq_name, "down", "others")


##QC

DMSO_AUXIN_horizontal_filter <-filter(DMSO_AUXIN_horizontal, (Rate01MAD.x/Rate01Median.x) < 0.75,  (Rate01MAD.y/Rate01Median.y) < 0.75 ,  
                                      (Rate10MAD.x/Rate10Median.x) < 0.75,  (Rate10MAD.y/Rate10Median.y) < 0.75 , 
                                      (EjectMAD.x/EjectMedian.x) < 0.75,  (EjectMAD.y/EjectMedian.y) < 0.75, 
                                      (BurstMAD.x/BurstMedian.x) < 0.75, (BurstMAD.y/BurstMedian.y) < 0.75,    Expression.x > 0.01,  Expression.y > 0.01)


##

DMSO_AUXIN_horizontal_filter <- transform(DMSO_AUXIN_horizontal_filter, LFC_off = log2((1/Rate01Median.y)/(1/Rate01Median.x)))

DMSO_AUXIN_horizontal_filter <- transform(DMSO_AUXIN_horizontal_filter, LFC_BS = log2(BurstMedian.y/BurstMedian.x))


##

DMSO_AUXIN_horizontal_filter <- transform(DMSO_AUXIN_horizontal_filter, OFF_time_DMSO = (1/Rate01Median.x))

DMSO_AUXIN_horizontal_filter <- transform(DMSO_AUXIN_horizontal_filter, OFF_time_MED26 = (1/Rate01Median.y))

DMSO_AUXIN_horizontal_filter <- transform(DMSO_AUXIN_horizontal_filter, BS_DMSO = BurstMedian.x)

DMSO_AUXIN_horizontal_filter <- transform(DMSO_AUXIN_horizontal_filter, BS_MED26 = BurstMedian.y)



##

ggplot(DMSO_AUXIN_horizontal_filter, aes(x=cat_IAA_120_bulk, y= log2(Expression.y/Expression.x), color=cat_IAA_120_bulk))  + 
  geom_boxplot(width = 0.3) + theme_classic2()  + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +   
  scale_color_brewer(palette="Set1") +  theme(axis.text=element_text(size=12),
                                              axis.title=element_text(size=14,face="bold"))  + geom_hline(yintercept=0, size = 1, linetype = 'dotted', col = 'black') +
  ylab("log2 (MED26/DMSO) Expression") + ylim(-4, 5) 


ggplot(DMSO_AUXIN_horizontal_filter, aes(x=cat_IAA_120_bulk, y= LFC_off, color=cat_IAA_120_bulk))  +   geom_boxplot(width = 0.3) + theme_classic2()  + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +   
  scale_color_brewer(palette="Set1") +  theme(axis.text=element_text(size=12),
                                              axis.title=element_text(size=14,face="bold"))  + geom_hline(yintercept=0, size = 1, linetype = 'dotted', col = 'black') + 
  ylab("log2 (MED26/DMSO) OFF time") + ylim(-4, 4)


ggplot(DMSO_AUXIN_horizontal_filter, aes(x=cat_IAA_120_bulk, y= LFC_BS, color=cat_IAA_120_bulk))   + 
  geom_boxplot(width = 0.3) + theme_classic2()  + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +   
  scale_color_brewer(palette="Set1") +  theme(axis.text=element_text(size=12),
                                              axis.title=element_text(size=14,face="bold"))  + geom_hline(yintercept=0, size = 1, linetype = 'dotted', col = 'black') + 
  ylab("log2 (MED26/DMSO) Burst Size") + ylim(-4, 4)
##

table(DMSO_AUXIN_horizontal_filter$cat_IAA_120_bulk)

down_median <- filter(DMSO_AUXIN_horizontal_filter, cat_IAA_120_bulk == "down")

median(down_median$LFC_off)

median(down_median$LFC_BS)

wilcox.test(down_median$OFF_time_DMSO , down_median$OFF_time_MED26, paired=TRUE)

wilcox.test(down_median$BS_DMSO , down_median$BS_MED26, paired=TRUE)

