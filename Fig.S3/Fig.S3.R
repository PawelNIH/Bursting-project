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
G1 <- read.csv(file = 'r1_G1_Summary_scRNA-ss-HCT116ns_2.csv')

G1_burst <- read.csv(file = 'r1_G1burst_scRNA-ss-HCT116ns_HCTDMSOr1G1_2.csv')


G1 <- merge(G1, G1_burst, by="Gene") 

G1$Cell_cycle <- "G1"

##

S <- read.csv(file = 'r1_S_Summary_scRNA-ss-HCT116ns_2.csv')

S_burst <- read.csv(file = 'r1_S_burst_scRNA-ss-HCT116ns_HCTDMSOr1S_2.csv')


S <- merge(S, S_burst, by="Gene")

S$Cell_cycle <- "S"

##
horizontal <- merge(G1, S, by="Gene") 

##

G2M <- read.csv(file = 'r1_G2M_Summary_scRNA-ss-HCT116ns_2.csv')

G2M_burst <- read.csv(file = 'r1_G2M_burst_scRNA-ss-HCT116ns_HCTDMSOr1G2M_2.csv')


G2M <- merge(G2M, G2M_burst, by="Gene") 

G2M$Cell_cycle <- "G2M"

##
horizontal <- merge(horizontal, G2M, by="Gene") 



horizontal_filtered <-filter(horizontal, (Rate01MAD.x/Rate01Median.x) < 0.75,  (Rate01MAD.y/Rate01Median.y) < 0.75,  (Rate01MAD/Rate01Median) < 0.75,
                            (Rate10MAD.x/Rate10Median.x) < 0.75,  (Rate10MAD.y/Rate10Median.y) < 0.75,  (Rate10MAD/Rate10Median) < 0.75, 
                            (EjectMAD.x/EjectMedian.x) < 0.75,  (EjectMAD.y/EjectMedian.y) < 0.75,   (EjectMAD/EjectMedian) < 0.75,
                            (BurstMAD.x/BurstMedian.x) < 0.75, (BurstMAD.y/BurstMedian.y) < 0.75,  (BurstMAD/BurstMedian) < 0.75,
                            Expression.x > 0.01, Expression.y > 0.01,  Expression > 0.01)




## # S vs. G1 filtered

ggplot(horizontal_filtered, aes(x= log10(1/Rate01Median.x), y= log10(1/Rate01Median.y))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "pearson", size = 5) 



ggplot(horizontal_filtered, aes(x= log10(BurstMedian.x), y= log10(BurstMedian.y))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "pearson", size = 5) 



# G2 vs. G1 filtered


ggplot(horizontal_filtered, aes(x= log10(1/Rate01Median.x), y= log10(1/Rate01Median))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "pearson", size = 5) 


ggplot(horizontal_filtered, aes(x= log10(BurstMedian.x), y= log10(BurstMedian))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "pearson", size = 5) 

