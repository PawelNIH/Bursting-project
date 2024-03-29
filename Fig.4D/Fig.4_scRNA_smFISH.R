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

SC_FISH<- read.csv(file = 'Fig.4_scRNA_smFISH.csv')


# #QC

SC_FISH <-filter(SC_FISH, (Rate01MAD/Rate01Median) < 0.75,  
            (Rate10MAD/Rate10Median) < 0.75,  
            (EjectMAD/EjectMedian) < 0.75,  
            (BurstMAD/BurstMedian)< 0.75,  
            Expression > 0.01)

## Add OFF time
SC_FISH <- transform(SC_FISH, off_time = 1/Rate01Median)

## Add error for OFF time
SC_FISH <- transform(SC_FISH, off_time_MAD = (1/Rate01Median^2*Rate01MAD))


## Add expression SD
SC_FISH <- transform(SC_FISH, expr_sd = sqrt(Variance))


### Burst size analysis
fishdmso_burst <- SC_FISH$BurstMedian[SC_FISH$Method=="smRNA_FISH" & SC_FISH$Condition=="HCTDMSO2h"]

fishauxin_burst <-SC_FISH$BurstMedian[SC_FISH$Method=="smRNA_FISH" & SC_FISH$Condition=="HCTAUXIN2h"]

logfish_burst = log2(fishauxin_burst) - log2(fishdmso_burst)

##

fishdmso_burstmad <- SC_FISH$BurstMAD[SC_FISH$Method=="smRNA_FISH" & SC_FISH$Condition=="HCTDMSO2h"]

fishauxin_burstmad <-SC_FISH$BurstMAD[SC_FISH$Method=="smRNA_FISH" & SC_FISH$Condition=="HCTAUXIN2h"]

fish_burst_error <- fishauxin_burstmad/fishauxin_burst + fishdmso_burstmad/fishdmso_burst


##
scrnadmso_burst <-SC_FISH$BurstMedian[SC_FISH$Method=="scRNA_Seq" & SC_FISH$Condition=="HCTDMSO2h"]

scrnaauxin_burst <-SC_FISH$BurstMedian[SC_FISH$Method=="scRNA_Seq" & SC_FISH$Condition=="HCTAUXIN2h"]

logscrna_burst = log2(scrnaauxin_burst) - log2(scrnadmso_burst)


##
scrnadmso_burstmad <-SC_FISH$BurstMAD[SC_FISH$Method=="scRNA_Seq" & SC_FISH$Condition=="HCTDMSO2h"]

scrnaauxin_burstmad <-SC_FISH$BurstMAD[SC_FISH$Method=="scRNA_Seq" & SC_FISH$Condition=="HCTAUXIN2h"]

scrna_burst_error = scrnaauxin_burstmad/scrnaauxin_burst + scrnadmso_burstmad/scrnadmso_burst


## Test model

model <- lm(formula = logscrna_burst ~ logfish_burst)

summary(model)


## Create a scatter plot with error bars 
plot(logfish_burst, logscrna_burst, ylim = c(-2, 2))
arrows(
  logfish_burst - fish_burst_error, logscrna_burst,
  logfish_burst + fish_burst_error, logscrna_burst,
  angle = 90, code = 0, length = 0.1, col = "black"
)
arrows(
  logfish_burst, logscrna_burst - scrna_burst_error,
  logfish_burst, logscrna_burst + scrna_burst_error,
  angle = 90, code = 0, length = 0.1, col = "black"
)

## Add ablines

abline(h = 0, col = "black", lwd = 1, lty = 2)
abline(v = 0, col = "black", lwd = 1, lty = 2)


## OFF time analysis



fishdmso_off <- SC_FISH$off_time[SC_FISH$Method=="smRNA_FISH" & SC_FISH$Condition=="HCTDMSO2h"]

fishauxin_off <-SC_FISH$off_time[SC_FISH$Method=="smRNA_FISH" & SC_FISH$Condition=="HCTAUXIN2h"]

logfish_off = log2(fishauxin_off) - log2(fishdmso_off)

##

fishdmso_offmad <- SC_FISH$off_time_MAD[SC_FISH$Method=="smRNA_FISH" & SC_FISH$Condition=="HCTDMSO2h"]

fishauxin_offmad <-SC_FISH$off_time_MAD[SC_FISH$Method=="smRNA_FISH" & SC_FISH$Condition=="HCTAUXIN2h"]

fish_off_error <- fishauxin_offmad/fishauxin_off + fishdmso_offmad/fishdmso_off


##
scrnadmso_off <-SC_FISH$off_time[SC_FISH$Method=="scRNA_Seq" & SC_FISH$Condition=="HCTDMSO2h"]

scrnaauxin_off <-SC_FISH$off_time[SC_FISH$Method=="scRNA_Seq" & SC_FISH$Condition=="HCTAUXIN2h"]

logscrna_off = log2(scrnaauxin_off) - log2(scrnadmso_off)


##
scrnadmso_offmad <-SC_FISH$off_time_MAD[SC_FISH$Method=="scRNA_Seq" & SC_FISH$Condition=="HCTDMSO2h"]

scrnaauxin_offmad <-SC_FISH$off_time_MAD[SC_FISH$Method=="scRNA_Seq" & SC_FISH$Condition=="HCTAUXIN2h"]

scrna_off_error = scrnaauxin_offmad/scrnaauxin_off + scrnadmso_offmad/scrnadmso_off



## Test model
model <- lm(formula = logscrna_off ~ logfish_off)

summary(model)



## Create a scatter plot with error bars 
plot(logfish_off, logscrna_off, xlim = c(-0.2, 0.7), ylim = c(-0.6, 1.6))

## Add ablines
abline(h = 0, col = "black", lwd = 1, lty = 2)
abline(v = 0, col = "black", lwd = 1, lty = 2)

## Add arrows 
arrows(
  logfish_off - fish_off_error, logscrna_off,
  logfish_off + fish_off_error, logscrna_off,
  angle = 90, code = 0, length = 0.1, col = "black"
)
arrows(
  logfish_off, logscrna_off - scrna_off_error,
  logfish_off, logscrna_off + scrna_off_error,
  angle = 90, code = 0, length = 0.1, col = "black"
)


## Expression


fishdmso_expr <- SC_FISH$Expression[SC_FISH$Method=="smRNA_FISH" & SC_FISH$Condition=="HCTDMSO2h"]

fishauxin_expr <-SC_FISH$Expression[SC_FISH$Method=="smRNA_FISH" & SC_FISH$Condition=="HCTAUXIN2h"]

logfish_expr = log2(fishauxin_expr) - log2(fishdmso_expr)

#

fishdmso_exprsd <- SC_FISH$expr_sd[SC_FISH$Method=="smRNA_FISH" & SC_FISH$Condition=="HCTDMSO2h"]

fishdmso_Gene <- SC_FISH$Gene[SC_FISH$Method=="smRNA_FISH" & SC_FISH$Condition=="HCTDMSO2h"]

fishdmso_no_cells <- c(3595, 4759, 4759, 4169, 5578, 2732, 5578)

fishdmso_exprsem <- fishdmso_exprsd / sqrt(fishdmso_no_cells)


##
fishauxin_exprsd <-SC_FISH$expr_sd[SC_FISH$Method=="smRNA_FISH" & SC_FISH$Condition=="HCTAUXIN2h"]

fishauxin_Gene <-SC_FISH$Gene[SC_FISH$Method=="smRNA_FISH" & SC_FISH$Condition=="HCTAUXIN2h"]

fishauxin_no_cells <- c(4523, 5382, 5382, 4523, 4591, 3453, 3921)

fishauxin_exprsem <- fishauxin_exprsd / sqrt(fishauxin_no_cells)

##

fish_expr_error <- fishauxin_exprsem/fishauxin_expr + fishdmso_exprsem/fishdmso_expr


##
scrnadmso_expr <-SC_FISH$Expression[SC_FISH$Method=="scRNA_Seq" & SC_FISH$Condition=="HCTDMSO2h"]

scrnaauxin_expr <-SC_FISH$Expression[SC_FISH$Method=="scRNA_Seq" & SC_FISH$Condition=="HCTAUXIN2h"]

logscrna_expr = log2(scrnaauxin_expr) - log2(scrnadmso_expr)


##
scrnadmso_exprsd <-SC_FISH$expr_sd[SC_FISH$Method=="scRNA_Seq" & SC_FISH$Condition=="HCTDMSO2h"]

scrnadmso_exprsem <- scrnadmso_exprsd / sqrt(3731)



scrnaauxin_exprsd <-SC_FISH$expr_sd[SC_FISH$Method=="scRNA_Seq" & SC_FISH$Condition=="HCTAUXIN2h"]

scrnaauxin_exprsem <- scrnaauxin_exprsd / sqrt(3481)

##

scrna_expr_error = scrnaauxin_exprsem/scrnaauxin_expr + scrnadmso_exprsem/scrnadmso_expr


## Test model
model <- lm(formula = logscrna_expr ~ logfish_expr)

summary(model)




## Create a scatter plot with error bars
plot(logfish_expr, logscrna_expr)

## Add ablines
abline(h = 0, col = "black", lwd = 1, lty = 2)
abline(v = 0, col = "black", lwd = 1, lty = 2)

## Add arrows
arrows(
  logfish_expr - fish_expr_error, logscrna_expr,
  logfish_expr + fish_expr_error, logscrna_expr,
  angle = 90, code = 0, length = 0.1, col = "black"
)
arrows(
  logfish_expr, logscrna_expr - scrna_expr_error,
  logfish_expr, logscrna_expr + scrna_expr_error,
  angle = 90, code = 0, length = 0.1, col = "black"
)


## OFF time without two outliers: SOX9 and MYC

SC_FISH <- read.csv(file = 'Fig.4_scRNA_smFISH.csv')

## Add OFF time
SC_FISH <- transform(SC_FISH, off_time = 1/Rate01Median)

## Add error for OFF time
SC_FISH <- transform(SC_FISH, off_time_MAD = (1/Rate01Median^2*Rate01MAD))

## Filter genes

SC_FISH_2 <- filter(SC_FISH, Gene == "ERRFI1" | Gene == "HMMR" | Gene == "KIF2C" | Gene == "PPM1G" | Gene == "PSMD14")


## Analysis

fishdmso_off2 <- SC_FISH_2$off_time[SC_FISH_2$Method=="smRNA_FISH" & SC_FISH_2$Condition=="HCTDMSO2h"]

fishauxin_off2 <-SC_FISH_2$off_time[SC_FISH_2$Method=="smRNA_FISH" & SC_FISH_2$Condition=="HCTAUXIN2h"]

logfish_off2 = log2(fishauxin_off2) - log2(fishdmso_off2)

##

fishdmso_offmad2 <- SC_FISH_2$off_time_MAD[SC_FISH_2$Method=="smRNA_FISH" & SC_FISH_2$Condition=="HCTDMSO2h"]

fishauxin_offmad2 <-SC_FISH_2$off_time_MAD[SC_FISH_2$Method=="smRNA_FISH" & SC_FISH_2$Condition=="HCTAUXIN2h"]

fish_off_error2 <- fishauxin_offmad2/fishauxin_off2 + fishdmso_offmad2/fishdmso_off2


##
scrnadmso_off2 <-SC_FISH_2$off_time[SC_FISH_2$Method=="scRNA_Seq" & SC_FISH_2$Condition=="HCTDMSO2h"]

scrnaauxin_off2 <-SC_FISH_2$off_time[SC_FISH_2$Method=="scRNA_Seq" & SC_FISH_2$Condition=="HCTAUXIN2h"]

logscrna_off2 = log2(scrnaauxin_off2) - log2(scrnadmso_off2)


##
scrnadmso_offmad2 <-SC_FISH_2$off_time_MAD[SC_FISH_2$Method=="scRNA_Seq" & SC_FISH_2$Condition=="HCTDMSO2h"]

scrnaauxin_offmad2 <-SC_FISH_2$off_time_MAD[SC_FISH_2$Method=="scRNA_Seq" & SC_FISH_2$Condition=="HCTAUXIN2h"]

scrna_off_error2 = scrnaauxin_offmad2/scrnaauxin_off2 + scrnadmso_offmad2/scrnadmso_off2


## Test model
model <- lm(formula = logscrna_off2 ~ logfish_off2)

summary(model)


## Create a scatter plot with error bars
plot(logfish_off2, logscrna_off2, xlim = c(-0.2, 0.7), ylim = c(-0.6, 1.2))

## Add ablines

abline(h = 0, col = "black", lwd = 1, lty = 2)
abline(v = 0, col = "black", lwd = 1, lty = 2)

## Add arrows 
arrows(
  logfish_off2 - fish_off_error2, logscrna_off2,
  logfish_off2 + fish_off_error2, logscrna_off2,
  angle = 90, code = 0, length = 0.1, col = "black"
)
arrows(
  logfish_off2, logscrna_off2 - scrna_off_error2,
  logfish_off2, logscrna_off2 + scrna_off_error2,
  angle = 90, code = 0, length = 0.1, col = "black"
)

