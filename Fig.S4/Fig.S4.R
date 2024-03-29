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

##Load data from Johnsson et al., Transcriptional kinetics and molecular functions of long noncoding RNAs, Nature Genetics volume 54, pages 306–317 (2022)

Nature <-read.csv(file = 'NatGen.csv')

## Filter genes that passed QC
Nature <- filter(Nature, pass.size == "TRUE", pass.kon == "TRUE", pass.mean == "TRUE")

##Load StochasticGene rates 
Fibro_CC <- read.csv(file = 'Summary_scRNA-ss-fibs_2.csv')

##Load StochasticGene burst size

Fibro_CC_burst <- read.csv("burst_scRNA-ss-fibs_CAST_2.csv")

##Merge data
Fibro_CC <- merge(Fibro_CC, Fibro_CC_burst, by="Gene")

##QC filtering (this paper)

Fibro_CC_filtered <- filter(Fibro_CC, (Rate01MAD/Rate01Median) < 0.75, 
                    (Rate10MAD/Rate10Median) < 0.75,  
                    (EjectMAD/EjectMedian) < 0.75,  
                    (BurstMAD/BurstMedian) < 0.75, 
                    Expression > 0.01)


##Comnpare StochasticGene and Johnsson et al. rates

Fibro_CC_filtered <- select(Fibro_CC_filtered, Gene, Rate01, Eject, Rate10, Rate01Median, EjectMedian, Rate10Median, Decay, Expression, BurstMedian)

colnames(Fibro_CC_filtered) <- c("GeneStableID", "Rate01", "Eject", "Rate10", "Rate01Median", "EjectMedian", "Rate10Median", "Decay", "Expression", "BurstMedian")


Nat_CC <- merge(Nature, Fibro_CC_filtered, by="GeneStableID") 

ggplot(Nat_CC, aes(x= log10(mean), y= log10(Expression))) + 
  theme_classic() +
  geom_point(size = 1
  ) + stat_cor(method = "spearman") + xlab("Expression Johnsson et al. (log10)") + ylab("Expression StochasticGene (log10)") +
  theme( axis.title = element_text(size = 7),  
    axis.text.x = element_text(size = 6),  
    axis.text.y = element_text(size = 6))

##StochasticGene posterior median rates

plot1 <- ggplot(Nat_CC, aes(x= log10(kon), y= log10(Rate01Median/Decay))) + 
  theme_classic() +
  geom_point(size = 1,  color = "gray"
  ) + stat_cor(method = "spearman") + xlab("ON rate Johnsson et al. (log10)") + ylab("ON rate Median (log10)")+
  theme( axis.title = element_text(size = 7),  
         axis.text.x = element_text(size = 6),  
         axis.text.y = element_text(size = 6))


plot2 <- ggplot(Nat_CC, aes(x= log10(koff), y= log10(Rate10Median/Decay))) + 
  theme_classic() +
  geom_point(size = 1,  color = "gray"
  ) + stat_cor(method = "spearman") + xlab("OFF rate Johnsson et al. (log10)") + ylab("OFF rate Median (log10)")+
  theme( axis.title = element_text(size = 7),  
         axis.text.x = element_text(size = 6),  
         axis.text.y = element_text(size = 6))

plot3 <- ggplot(Nat_CC, aes(x= log10(ksyn), y= log10(EjectMedian/Decay))) + 
  theme_classic() +
  geom_point(size = 1,  color = "gray"
  ) + stat_cor(method = "spearman") + xlab("Eject Johnsson et al. (log10)") + ylab("Eject Median (log10)")+
  theme( axis.title = element_text(size = 7),  
         axis.text.x = element_text(size = 6),  
         axis.text.y = element_text(size = 6))


plot4 <- ggplot(Nat_CC, aes(x= log10(ksyn/koff), y= log10(BurstMedian/Decay))) + 
  theme_classic() +
  geom_point(size = 1,  color = "gray"
  ) + stat_cor(method = "spearman") + xlab("Burst Size Johnsson et al. (log10)") + ylab("Burst Size Median (log10)")+
  theme( axis.title = element_text(size = 7),  
         axis.text.x = element_text(size = 6),  
         axis.text.y = element_text(size = 6))

ggarrange(plot1, plot2, plot3, plot4, ncol = 2, nrow = 2)



##StochasticGene MLE rates
plot5 <- ggplot(Nat_CC, aes(x= log10(kon), y= log10(Rate01/Decay))) + 
  theme(axis.title = element_text(size = 20)) +  theme(axis.text.x = element_text(size = 20)) +     theme(axis.text.y = element_text(size = 20))  +
  theme_classic() +
  geom_point(size = 1,  color = "gray"
  ) + stat_cor(method = "spearman") + xlab("ON rate Johnsson et al. (log10)") + ylab("ON rate MLE (log10)")

plot6 <- ggplot(Nat_CC, aes(x= log10(koff), y= log10(Rate10/Decay))) + 
  theme_classic() +
  geom_point(size = 1,  color = "gray"
  ) + stat_cor(method = "spearman") + xlab("OFF rate Johnsson et al. (log10)") + ylab("OFF rate MLE (log10)")

plot7 <- ggplot(Nat_CC, aes(x= log10(ksyn), y= log10(Eject/Decay))) + 
  theme_classic() +
  geom_point(size = 1,  color = "gray"
  ) + stat_cor(method = "spearman")+ xlab("Eject Johnsson et al. (log10)") + ylab("Eject MLE (log10)")

plot8 <- ggplot(Nat_CC, aes(x= log10(ksyn/koff), y= log10((Eject/Rate10))/Decay)) + 
  theme_classic() +
  geom_point(size = 1,  color = "gray"
  ) + stat_cor(method = "spearman") + xlab("Burst Size Johnsson et al. (log10)") + ylab("Burst Size MLE (log10)")

ggarrange(plot5, plot6, plot7, plot8, ncol = 2, nrow = 2)

### Plot histograms MAD/rate for genes that passed QC

histo_filtered <- filter(Fibro_CC, (Rate01MAD/Rate01Median) < 0.75, 
                      (Rate10MAD/Rate10Median) < 0.75,  
                      (EjectMAD/EjectMedian) < 0.75,  
                      (BurstMAD/BurstMedian) < 0.75, 
                      Expression > 0.01)

xlim <- c(0, 0.8)
num_bins <- 50
par(mfrow = c(1, 4))

a <- hist((histo_filtered$Rate01MAD / histo_filtered$Decay) / (histo_filtered$Rate01Median / histo_filtered$Decay), 
          col = "blue", 
          border = "blue", 
          xlim = xlim,
          breaks = seq(xlim[1], xlim[2], length.out = num_bins + 1), 
          main = "ON MAD/
          ON rate",   xlab = "")

b <- hist((histo_filtered$Rate10MAD / histo_filtered$Decay) / (histo_filtered$Rate10Median / histo_filtered$Decay), 
          col = "blue", 
          border = "blue", 
          xlim = xlim,
          breaks = seq(xlim[1], xlim[2], length.out = num_bins + 1), 
          main = "OFF MAD/
          OFF rate",  xlab = "")

c <- hist((histo_filtered$EjectMAD / histo_filtered$Decay) / (histo_filtered$EjectMedian / histo_filtered$Decay), 
          col = "blue", 
          border = "blue", 
          xlim = xlim,
          breaks = seq(xlim[1], xlim[2], length.out = num_bins + 1), 
          main = "Eject MAD/
          Eject rate",  xlab = "")

d <- hist((histo_filtered$BurstMAD / histo_filtered$Decay) / (histo_filtered$BurstMedian / histo_filtered$Decay), 
          col = "blue", 
          border = "blue", 
          xlim = xlim,
          breaks = seq(xlim[1], xlim[2], length.out = num_bins + 1), 
          main = "Burst Size MAD/
          Burst Size",  xlab = "")
 


##Plot correlations between rates
##Load correlation values  StochasticGene for data from Johnsson et al.

CAST_CC <- read.csv(file = "CAST_CC.csv", header = F)

colnames(CAST_CC) <- c("Gene", "x", "ON_OFF", "ON_eject", "OFF_eject")

##Load genes that passed QC
CAST_QC <- read.csv(file = "Fibro_QC.csv")

CAST_QC <- select(CAST_QC, Gene)

##Filter only genes that passed QC
CAST_CC_QC <- merge(CAST_CC, CAST_QC, by="Gene") 

##Plot
num_bins <- 50

par(mfrow = c(1, 3))

#plot the first histogram 


hist(CAST_CC_QC$ON_OFF, 
     col = "blue", 
     border = "blue",  
     breaks = seq(min(CAST_CC_QC$ON_OFF), max(CAST_CC_QC$ON_OFF), length.out = num_bins + 1), 
     main = "ON-OFF rate
     corr",  xlab = "")

xlim(c(0, 1))


#plot the second histogram 

hist(CAST_CC_QC$ON_eject, 
     col = "blue", 
     border = "blue",  
     breaks = seq(min(CAST_CC_QC$ON_eject), max(CAST_CC_QC$ON_eject), length.out = num_bins + 1), 
     main = "ON-Eject rate
     corr",  xlab = "")

xlim(c(0, 1))


# plot the third histogram 

hist(CAST_CC_QC$OFF_eject, 
     col = "blue", 
     border = "blue",  
     breaks = seq(min(CAST_CC_QC$OFF_eject), max(CAST_CC_QC$OFF_eject), length.out = num_bins + 1),
     main = "OFF-Eject rate
     corr",  xlab = "")
xlim(c(0, 1))


