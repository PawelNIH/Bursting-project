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

##Load summary data
Data <- read.csv(file = 'X')

##Load data for burst size estimations

Data_burst <- read.csv(file = 'Y')

##Merge data
Data <- merge(Data, Data_burst, by="Gene")

## QC filter


Data_filtered <- filter(Data, (Rate01MAD/Rate01Median) < 0.75, 
                    (Rate10MAD/Rate10Median) < 0.75,  
                    (EjectMAD/EjectMedian) < 0.75,  
                    (BurstMAD/BurstMedian) < 0.75, 
                    Expression > 0.01)
