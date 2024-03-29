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


##Load data WT
WT <- read.csv(file = 'Summary_scRNA-ss-aBn_2_WTr1.csv')

WT_burst <- read.csv(file = 'burst_scRNA-ss-aBn_xWTr1_2.csv')

WT <- merge(WT, WT_burst, by="Gene") 



##Load data KO

MYC <- read.csv(file = 'Summary_scRNA-ss-aBn_2_Mycr1.csv')

MYC_burst <- read.csv(file = 'burst_scRNA-ss-aBn_Mycmr1_2.csv')

MYC <- merge(MYC, MYC_burst, by="Gene") 


##Merge data
DMSO_AUXIN_horizontal <-  merge(WT, MYC, by="Gene") 


##QC filtering
DMSO_AUXIN_horizontal_filter <-filter(DMSO_AUXIN_horizontal, (Rate01MAD.x/Rate01Median.x) < 0.75,  (Rate01MAD.y/Rate01Median.y) < 0.75 ,  
                                      (Rate10MAD.x/Rate10Median.x) < 0.75,  (Rate10MAD.y/Rate10Median.y) < 0.75 , 
                                      (EjectMAD.x/EjectMedian.x) < 0.75,  (EjectMAD.y/EjectMedian.y) < 0.75, 
                                      (BurstMAD.x/BurstMedian.x) < 0.75, (BurstMAD.y/BurstMedian.y) < 0.75,    Expression.x > 0.01,  Expression.y > 0.01)



##DEG based on scRNA-Seq

DEG <- read.table("Mycmm_DE2_DEseq.txt", header=T)

down_DEG <- filter(DEG, avg_log2FC < -0.32)


DMSO_AUXIN_horizontal_filter$cat_DEG <- ifelse(DMSO_AUXIN_horizontal_filter$Gene %in% down_DEG$Gene, "down", "other")


DMSO_AUXIN_horizontal_filter <- transform(DMSO_AUXIN_horizontal_filter, LFC_off = log2((1/Rate01Median.y)/(1/Rate01Median.x)))

DMSO_AUXIN_horizontal_filter <- transform(DMSO_AUXIN_horizontal_filter, LFC_BS = log2(BurstMedian.y/BurstMedian.x))

DMSO_AUXIN_horizontal_filter <- transform(DMSO_AUXIN_horizontal_filter, OFF_time_WT = (1/Rate01Median.x))

DMSO_AUXIN_horizontal_filter <- transform(DMSO_AUXIN_horizontal_filter, OFF_time_KO = (1/Rate01Median.y))

DMSO_AUXIN_horizontal_filter <- transform(DMSO_AUXIN_horizontal_filter, BS_WT = BurstMedian.x)

DMSO_AUXIN_horizontal_filter <- transform(DMSO_AUXIN_horizontal_filter, BS_KO = BurstMedian.y)


##Plots
DMSO_AUXIN_horizontal_filter <- transform(DMSO_AUXIN_horizontal_filter, BS = log2(BurstMedian.y/BurstMedian.x))

 ggplot(DMSO_AUXIN_horizontal_filter, aes(x=cat_DEG, y= log2(Expression.y/Expression.x), color=cat_DEG))  + 
  geom_boxplot(width = 0.3) + theme_classic2()  + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +   
  scale_color_brewer(palette="Set1") +  theme(axis.text=element_text(size=12),
                                              axis.title=element_text(size=14,face="bold"))  + geom_hline(yintercept=0, size = 1, linetype = 'dotted', col = 'black') + ylim(-4, 4)


ggplot(DMSO_AUXIN_horizontal_filter, aes(x=cat_DEG, y= LFC_off, color=cat_DEG))  + 
  geom_boxplot(width = 0.3) + theme_classic2()  + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +   
  scale_color_brewer(palette="Set1") +  theme(axis.text=element_text(size=12),
                                              axis.title=element_text(size=14,face="bold"))  + geom_hline(yintercept=0, size = 1, linetype = 'dotted', col = 'black') + ylim(-4, 4) +
  ylab("log2 (KO/WT) OFF time") 


ggplot(DMSO_AUXIN_horizontal_filter, aes(x=cat_DEG, y= LFC_BS, color=cat_DEG))   + 
  geom_boxplot(width = 0.3) + theme_classic2()  + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +   
  scale_color_brewer(palette="Set1") +  theme(axis.text=element_text(size=12),
                                              axis.title=element_text(size=14,face="bold"))  + geom_hline(yintercept=0, size = 1, linetype = 'dotted', col = 'black') + ylim(-4, 4) +
  ylab("log2 (KO/WT) Burst Size") 



###

table(DMSO_AUXIN_horizontal_filter$cat_DEG)

down_median <- filter(DMSO_AUXIN_horizontal_filter, cat_DEG == "down")

median(down_median$LFC_off)

median(down_median$LFC_BS)

wilcox.test(down_median$OFF_time_WT , down_median$OFF_time_KO, paired=TRUE)

wilcox.test(down_median$BS_WT , down_median$BS_KO, paired=TRUE)

