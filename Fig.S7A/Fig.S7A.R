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



## RAD21 degron


RAD21_bulk_RNA <- read.table("HCT116_RAD21_degron_fdroutput_list", sep=",", header=T)

RAD21_bulk_RNA <- filter(RAD21_bulk_RNA, gt=="mRNAseq_120m")

down_120 <- subset(RAD21_bulk_RNA, gt=="mRNAseq_120m" & ud == "down")
up_120 <- subset(RAD21_bulk_RNA, gt=="mRNAseq_120m" & ud == "up")

RAD21_bulk_RNA$cat_IAA_120_bulk <- ifelse(RAD21_bulk_RNA$refseq_name %in% down_120$refseq_name, "down", "others")



## JQ1


JQ1_bulk_RNA <- read.table("HCT_JQ1_fdroutput_list", sep=",", header=T)


down_120 <- subset(JQ1_bulk_RNA,  ud == "down")
up_120 <- subset(JQ1_bulk_RNA, ud == "up")

JQ1_bulk_RNA$cat_IAA_120_bulk <- ifelse(JQ1_bulk_RNA$refseq_name %in% down_120$refseq_name, "down", "others")



## MED26 degron
MED26_bulk_RNA <- read.table("HCT_Med26Degron_fdroutput_list", sep=",", header=T)


down_120 <- subset(MED26_bulk_RNA,  ud == "down")
up_120 <- subset(MED26_bulk_RNA, ud == "up")

MED26_bulk_RNA$cat_IAA_120_bulk <- ifelse(MED26_bulk_RNA$refseq_name %in% down_120$refseq_name, "down", "others")


## combine


RAD21_bulk_RNA$Treatment <- "-RAD21"

JQ1_bulk_RNA$Treatment <- "JQ1"

JQ1_bulk_RNA$gt <- "mRNAseq_120m"


MED26_bulk_RNA$Treatment <- "-MED26"

MED26_bulk_RNA$gt <- "mRNAseq_120m"


combined_dataframe <- rbind(RAD21_bulk_RNA, JQ1_bulk_RNA, MED26_bulk_RNA)

combined_dataframe <- filter(combined_dataframe, cat_IAA_120_bulk == "down")

combined_dataframe$Treatment <- factor(combined_dataframe$Treatment, 
                                       levels = c("-RAD21", "JQ1", "-MED26"))

ggplot(combined_dataframe, aes(x=Treatment, y= log2FoldChange, fill=Treatment)) + 
  geom_boxplot(width = 0.3, outlier.colour = "black", outlier.shape = 16, outlier.size = 1) + 
  theme_classic2() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +   
  scale_fill_manual(values = c("red", "lightblue", "magenta")) +  
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14, face="bold")) +
  geom_hline(yintercept=0, size = 1, linetype = 'dotted', col = 'black') +
  ylim(-5, 2)

table(combined_dataframe$Treatment)


