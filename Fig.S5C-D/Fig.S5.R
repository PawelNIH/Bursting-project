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

## Load data
ES <- read.csv(file = "Summary_scRNA-ss-ESCDMSO_2.csv")

ES_MAD <- read.csv(file = 'burst_scRNA-ss-ESC_ESCCDMSO_2.csv')

ES <- merge(ES, ES_MAD, by="Gene") 

ES <- transform(ES, norm_expr = Expression/18786.05*16445)

ES <- transform(ES, norm_BS = BurstMedian/18786.05*16445)


aB <- read.csv(file = 'Summary_scRNA-ss-aBn_2.csv')


aB_MAD <- read.csv(file = 'burst_scRNA-ss-aBn_CtrlSpaB_2.csv')

aB <- merge(aB, aB_MAD, by="Gene") 

aB <- transform(aB, norm_expr = Expression/11940.82*16445)

aB <- transform(aB, norm_BS = BurstMedian/11940.82*16445)


ES_aB <- merge(ES, aB, by="Gene") 

#Filter

ES_aB <-filter(ES_aB, (Rate01MAD.x/Rate01Median.x) < 0.75,  (Rate01MAD.y/Rate01Median.y) < 0.75 ,  
                     (Rate10MAD.x/Rate10Median.x) < 0.75,  (Rate10MAD.y/Rate10Median.y) < 0.75 , 
                     (EjectMAD.x/EjectMedian.x) < 0.75,  (EjectMAD.y/EjectMedian.y) < 0.75, 
               (BurstMAD.x/BurstMedian.x) < 0.75, (BurstMAD.y/BurstMedian.y) < 0.75,    Expression.x > 0.01,  Expression.y > 0.01)




##ChIA_PET_cell_type_specific_enhancers_ES_B from Kieffer-Kwon et al., Cell 2013: Interactome maps of mouse gene regulatory domains reveal basic principles of transcriptional regulation

enhancer <- read.csv(file = "ChiAPET_ES_B.csv")

ES_aB_enhancer<- merge(ES_aB, enhancer, by="Gene") 

##more or less than 2 or 4 enhancers (edit number of enhancers: BminES)

ES_aB_enhancer <- transform(ES_aB_enhancer, BminES = B_enhancer - ES_enhancer)

lost <- filter(ES_aB_enhancer, Type == "> ES cells" & BminES <= -4)

gained <- filter(ES_aB_enhancer, Type == "> B cells" & BminES >= 4)
 


ES_aB_enhancer$Test <- ifelse(ES_aB_enhancer$Gene %in% lost$Gene, "lost", "no")

ES_aB_enhancer$Test <- ifelse(ES_aB_enhancer$Gene %in% gained$Gene, "gained", ES_aB_enhancer$Test)
ES_aB_enhancer <- filter(ES_aB_enhancer, Test != "no")


##plots


ES_aB_enhancer <- transform(ES_aB_enhancer, LFC_exp = log2((norm_expr.y)/(norm_expr.x)))


ES_aB_enhancer <- transform(ES_aB_enhancer, LFC_off = log2((1/Rate01Median.y)/(1/Rate01Median.x)))


ES_aB_enhancer <- transform(ES_aB_enhancer, LFC_BS = log2(norm_BS.y/norm_BS.x))



## 
my_comparisons <- list(c("gained", "lost"))

ggplot(ES_aB_enhancer, aes(x=Test, y= LFC_exp, color=Test))+ 
  geom_boxplot(width = 0.2) + theme_classic2()  + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  scale_color_brewer(palette="Set1") +  theme(axis.text=element_text(size=12),
                                              axis.title=element_text(size=14,face="bold")) +
  stat_compare_means(comparisons=my_comparisons, method = "wilcox", aes(label=..p.adj..)) + ylim(-5, 7)

wilcox.test(ES_aB_enhancer$LFC_exp ~ ES_aB_enhancer$Test)


##
ggplot(ES_aB_enhancer, aes(x=Test, y= LFC_off, color=Test))+ 
  geom_boxplot(width = 0.2) + theme_classic2()  + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  scale_color_brewer(palette="Set1") +  theme(axis.text=element_text(size=12),
                                              axis.title=element_text(size=14,face="bold")) +
  stat_compare_means(comparisons=my_comparisons, method = "wilcox", aes(label=..p.adj..)) + ylim(-5, 7)


wilcox.test(ES_aB_enhancer$LFC_off ~ ES_aB_enhancer$Test)

##

ggplot(ES_aB_enhancer, aes(x=Test, y= LFC_BS, color=Test)) + 
  geom_boxplot(width = 0.2) + theme_classic2()  + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  labs( y = "log2 (B/ES BS)")+
  scale_color_brewer(palette="Set1") +  theme(axis.text=element_text(size=12),
                                              axis.title=element_text(size=14,face="bold")) +
  stat_compare_means(comparisons=my_comparisons, method = "wilcox", aes(label=..p.adj..))+ ylim(-5, 7)


wilcox.test(ES_aB_enhancer$LFC_BS ~ ES_aB_enhancer$Test)



