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

Pim1 <- read.csv(file = "Fig.3D-F_nonnormalized_rates.csv")

## Average reads per cell across analyzed cell types: 16445 (use for normalizatoon)



aB <- filter(Pim1, Cell == "aB") 

aB <- transform(aB, SD_expr = sqrt(Variance))

aB <- transform(aB, SEM_expr = SD_expr/sqrt(6453.00))

aB <- transform(aB, norm_eject = EjectMedian/11940.82*16445)

aB <- transform(aB, norm_expr = Expression/11940.82*16445)

aB <- transform(aB, norm_BS = BurstMedian/11940.82*16445)


aB <- transform(aB, norm_SEM_expr = SEM_expr/11940.82*16445)

aB <- transform(aB, norm_BS_MAD = BurstMAD/11940.82*16445)


##

ES <- filter(Pim1, Cell == "ESr1") 

ES <- transform(ES, SD_expr = sqrt(Variance))

ES <- transform(ES, SEM_expr = SD_expr/sqrt(6550.00))

ES <- transform(ES, norm_eject = EjectMedian/18786.05*16445)

ES <- transform(ES, norm_expr = Expression/18786.05*16445)

ES <- transform(ES, norm_BS = BurstMedian/18786.05*16445)


ES <- transform(ES, norm_SEM_expr = SEM_expr/18786.05*16445)

ES <- transform(ES, norm_BS_MAD = BurstMAD/18786.05*16445)


##
ESr2 <- filter(Pim1, Cell == "ESr2")

ESr2 <- transform(ESr2, SD_expr = sqrt(Variance))

ESr2 <- transform(ESr2, SEM_expr = SD_expr/sqrt(3593.00))

ESr2 <- transform(ESr2, norm_eject = EjectMedian/31926.75*16445)

ESr2 <- transform(ESr2, norm_expr = Expression/31926.75*16445)

ESr2 <- transform(ESr2, norm_BS = BurstMedian/31926.75*16445)


ESr2 <- transform(ESr2, norm_SEM_expr = SEM_expr/31926.75*16445)

ESr2 <- transform(ESr2, norm_BS_MAD = BurstMAD/31926.75*16445)


##

BMr1 <- filter(Pim1, Cell == "BMr1")

BMr1 <- transform(BMr1, SD_expr = sqrt(Variance))

BMr1 <- transform(BMr1, SEM_expr = SD_expr/sqrt(3969))

BMr1 <- transform(BMr1, norm_eject = EjectMedian/3713.008062*16445)

BMr1 <- transform(BMr1, norm_expr = Expression/3713.008062*16445)

BMr1 <- transform(BMr1, norm_BS = BurstMedian/3713.008062*16445)



BMr1 <- transform(BMr1, norm_SEM_expr = SEM_expr/3713.008062*16445)

BMr1 <- transform(BMr1, norm_BS_MAD = BurstMAD/3713.008062*16445)


##

BMr2 <- filter(Pim1, Cell ==  "BMr2")

BMr2 <- transform(BMr2, SD_expr = sqrt(Variance))

BMr2 <- transform(BMr2, SEM_expr = SD_expr/sqrt(1767))

BMr2 <- transform(BMr2, norm_eject = EjectMedian/4709.259196*16445)

BMr2 <- transform(BMr2, norm_expr = Expression/4709.259196*16445)

BMr2 <- transform(BMr2, norm_BS = BurstMedian/4709.259196*16445)


BMr2 <- transform(BMr2, norm_SEM_expr = SEM_expr/4709.259196*16445)

BMr2 <- transform(BMr2, norm_BS_MAD = BurstMAD/4709.259196*16445)


##

LNK <- filter(Pim1, Cell == "LNK")

LNK <- transform(LNK, SD_expr = sqrt(Variance))

LNK <- transform(LNK, SEM_expr = SD_expr/sqrt(4932))

LNK <- transform(LNK, norm_eject = EjectMedian/27295.22607*16445)

LNK <- transform(LNK, norm_expr = Expression/27295.22607*16445)

LNK <- transform(LNK, norm_BS = BurstMedian/27295.22607*16445)



LNK <- transform(LNK, norm_SEM_expr = SEM_expr/27295.22607*16445)

LNK <- transform(LNK, norm_BS_MAD = BurstMAD/27295.22607*16445)



##

MAST <- filter(Pim1, Cell == "MAST")

MAST <- transform(MAST, SD_expr = sqrt(Variance))

MAST <- transform(MAST, SEM_expr = SD_expr/sqrt(4279))

MAST <- transform(MAST, norm_eject = EjectMedian/6586.838747*16445)

MAST <- transform(MAST, norm_expr = Expression/6586.838747*16445)

MAST <- transform(MAST, norm_BS = BurstMedian/6586.838747*16445)



MAST <- transform(MAST, norm_SEM_expr = SEM_expr/6586.838747*16445)

MAST <- transform(MAST, norm_BS_MAD = BurstMAD/6586.838747*16445)



##

rB <- filter(Pim1, Cell == "rB")

rB <- transform(rB, SD_expr = sqrt(Variance))

rB <- transform(rB, SEM_expr = SD_expr/sqrt(1867.00))

rB <- transform(rB, norm_eject = EjectMedian/5571.08*16445)

rB <- transform(rB, norm_expr = Expression/5571.08*16445)

rB <- transform(rB, norm_BS = BurstMedian/5571.08*16445)




rB <- transform(rB, norm_SEM_expr = SEM_expr/5571.08*16445)

rB <- transform(rB, norm_BS_MAD = BurstMAD/5571.08*16445)



##

SKINr1 <- filter(Pim1, Cell == "SKINr1")

SKINr1 <- transform(SKINr1, SD_expr = sqrt(Variance))

SKINr1 <- transform(SKINr1, SEM_expr = SD_expr/sqrt(2306))

SKINr1 <- transform(SKINr1, norm_eject = EjectMedian/19426.4588*16445)

SKINr1 <- transform(SKINr1, norm_expr = Expression/19426.4588*16445)

SKINr1 <- transform(SKINr1, norm_BS = BurstMedian/19426.4588*16445)




SKINr1 <- transform(SKINr1, norm_SEM_expr = SEM_expr/19426.4588*16445)

SKINr1 <- transform(SKINr1, norm_BS_MAD = BurstMAD/19426.4588*16445)



##

SKINr2 <- filter(Pim1, Cell == "SKINr2")

SKINr2 <- transform(SKINr2, SD_expr = sqrt(Variance))

SKINr2 <- transform(SKINr2, SEM_expr = SD_expr/sqrt(1488))

SKINr2 <- transform(SKINr2, norm_eject = EjectMedian/20338.43212*16445)

SKINr2 <- transform(SKINr2, norm_expr = Expression/20338.43212*16445)

SKINr2 <- transform(SKINr2, norm_BS = BurstMedian/20338.43212*16445)



SKINr2 <- transform(SKINr2, norm_SEM_expr = SEM_expr/20338.43212*16445)

SKINr2 <- transform(SKINr2, norm_BS_MAD = BurstMAD/20338.43212*16445)


##

CH12 <- filter(Pim1, Cell  == "CH12")

CH12 <- transform(CH12, SD_expr = sqrt(Variance))

CH12 <- transform(CH12, SEM_expr = SD_expr/sqrt(7261))

CH12 <- transform(CH12, norm_eject = EjectMedian/7078.01873*16445)

CH12 <- transform(CH12, norm_expr = Expression/7078.01873*16445)

CH12 <- transform(CH12, norm_BS = BurstMedian/7078.01873*16445)




CH12 <- transform(CH12, norm_SEM_expr = SEM_expr/7078.01873*16445)

CH12 <- transform(CH12, norm_BS_MAD = BurstMAD/7078.01873*16445)



##

CH12r1 <- filter(Pim1, Cell ==  "CH12r1")

CH12r1 <- transform(CH12r1, SD_expr = sqrt(Variance))

CH12r1 <- transform(CH12r1, SEM_expr = SD_expr/sqrt(6273))

CH12r1 <- transform(CH12r1, norm_eject = EjectMedian/39963.34959*16445)

CH12r1 <- transform(CH12r1, norm_expr = Expression/39963.34959*16445)

CH12r1 <- transform(CH12r1, norm_BS = BurstMedian/39963.34959*16445)


CH12r1 <- transform(CH12r1, norm_SEM_expr = SEM_expr/39963.34959*16445)

CH12r1 <- transform(CH12r1, norm_BS_MAD = BurstMAD/39963.34959*16445)



###



all <- dplyr::bind_rows(ES, ESr2, aB, rB, MAST, SKINr1, SKINr2, LNK, CH12, CH12r1, BMr1, BMr2)



all_filter <-filter(all, (Rate01MAD/Rate01Median) < 0.75,  
                    (Rate10MAD/Rate10Median) < 0.75,  
                    (EjectMAD/EjectMedian) < 0.75,  
                    (BurstMAD/BurstMedian)< 0.75,  
                    Expression > 0.01)


Pim_filtered <- all_filter

Pim_filtered <- transform(Pim_filtered, off_time = (1/Rate01Median))

Pim_filtered <- transform(Pim_filtered, off_time_MAD = (1/Rate01Median^2*Rate01MAD))

#Fig.3D


##OFF time

OFF_time_aB <- Pim_filtered$off_time[Pim_filtered$Cell == "aB"]

OFF_timeMAD_aB <- Pim_filtered$off_time_MAD[Pim_filtered$Cell == "aB"]

##
OFF_time_ESr1 <- Pim_filtered$off_time[Pim_filtered$Cell == "ESr1"]

OFF_timeMAD_ESr1 <- Pim_filtered$off_time_MAD[Pim_filtered$Cell == "ESr1"]

##
OFF_time_ESr2 <- Pim_filtered$off_time[Pim_filtered$Cell == "ESr2"]

OFF_timeMAD_ESr2 <- Pim_filtered$off_time_MAD[Pim_filtered$Cell == "ESr2"]

##

OFF_time_ESr1_2 <- (OFF_time_ESr1+OFF_time_ESr2) / 2

##
FC_off_time <- log2(OFF_time_aB / OFF_time_ESr1_2)



##normalized burst size

BS_aB <- Pim_filtered$norm_BS[Pim_filtered$Cell == "aB"]

BS_MAD_aB <- Pim_filtered$norm_BS_MAD[Pim_filtered$Cell == "aB"]

##
BS_ESr1 <- Pim_filtered$norm_BS[Pim_filtered$Cell == "ESr1"]

BS_MAD_ESr1 <- Pim_filtered$norm_BS_MAD[Pim_filtered$Cell == "ESr1"]


##
BS_ESr2 <- Pim_filtered$norm_BS[Pim_filtered$Cell == "ESr2"]

BS_MAD_ESr2 <- Pim_filtered$norm_BS_MAD[Pim_filtered$Cell == "ESr2"]


##


BS_ESr1_2 <-(BS_ESr1+BS_ESr2) / 2

FC_BS <- log2(BS_aB / BS_ESr1_2)


##Fig. 3D

par(mfrow = c(2, 1))

###x

barplot(
  FC_off_time,
  beside = FALSE,
  col = c("gray", "lightblue"),
  main = "LFC OFF time",
  width = 0.01,
  ylim = c(-3, 3)
)

abline(h = 0, col = "black", lty = 2)

barplot(
  FC_BS,
  beside = FALSE,
  col = c("gray", "lightblue"),
  main = "LFC norm BS",
  width = 0.01,
  ylim = c(-3, 3)
)
abline(h = 0, col = "black", lty = 2)


#Fig. 3E-F

ggplot(Pim_filtered, aes(x= (norm_expr), y= (off_time)))  + 
  theme_classic() +
  geom_point(aes(colour = Cell), size = 5, 
  ) + 
  geom_smooth(method='lm', se = FALSE) +
  geom_errorbar(aes(xmin = norm_expr - norm_SEM_expr, xmax = norm_expr + norm_SEM_expr))+ 
  geom_errorbar(aes(ymin = off_time - off_time_MAD, ymax = off_time + off_time_MAD)) 


model <- lm(Pim_filtered$off_time~Pim_filtered$norm_expr)
summary(model)

ggplot(Pim_filtered, aes(x= (norm_expr), y= (norm_BS)))  + 
  theme_classic() +
  geom_point(aes(colour = Cell),size = 5, 
  ) + 
  geom_smooth(method='lm', se = FALSE)  +
  geom_errorbar(aes(xmin = norm_expr - norm_SEM_expr, xmax = norm_expr + norm_SEM_expr)) + 
  geom_errorbar(aes(ymin = norm_BS - norm_BS_MAD, ymax = norm_BS + norm_BS_MAD)) 

model2 <- lm(Pim_filtered$norm_BS~Pim_filtered$norm_expr)
summary(model2)

