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

##rep1 Load data (QC filtered rates)
rep1 <- read.csv(file = 'rep1_rates_HK_TF.csv')

rep1$rep <- "rep1"

##rep2 Load data (QC filtered rates)
rep2 <- read.csv(file = 'rep2_rates_HK_TF.csv')

rep2$rep <- "rep2"

##Combine replicates - plots without error bars

Compare <- dplyr::bind_rows(rep1, rep2)

a<- ggplot(Compare, aes(x = log10(Expression), color = Type, linetype = as.factor(rep))) + 
  stat_ecdf(geom = "step", size = 1) +   
  theme_classic2() + 
  scale_color_brewer(palette = "Set1") +
  scale_linetype_discrete(name = "Replicate") +
  theme(axis.text = element_text(size = 6), 
        axis.title = element_text(size = 7))

b<- ggplot(Compare, aes(x = log10(Decay), color = Type, linetype = as.factor(rep))) + 
  stat_ecdf(geom = "step", size = 1) +   
  theme_classic2() + 
  scale_color_brewer(palette = "Set1") +
  scale_linetype_discrete(name = "Replicate") +
  theme(axis.text = element_text(size = 6), 
        axis.title = element_text(size = 7))

c<- ggplot(Compare, aes(x = log10(1/Rate01Median), color = Type, linetype = as.factor(rep))) + 
  stat_ecdf(geom = "step", size = 1) +   
  theme_classic2() + 
  scale_color_brewer(palette = "Set1") +
  scale_linetype_discrete(name = "Replicate") +
  theme(axis.text = element_text(size = 6), 
        axis.title = element_text(size = 7))

d<- ggplot(Compare, aes(x = log10(BurstMedian), color = Type, linetype = as.factor(rep))) + 
  stat_ecdf(geom = "step", size = 1) +   
  theme_classic2() + 
  scale_color_brewer(palette = "Set1") +
  scale_linetype_discrete(name = "Replicate") +
  theme(axis.text = element_text(size = 6), 
        axis.title = element_text(size = 7)) +   theme_classic2() 

ggarrange(a,b, c, d, ncol = 2, nrow = 2)
table(Compare$rep, Compare$Type)


##Test Expression rep1
test_HK <- filter(Compare, Type == "HK" & rep == "rep1")

test_HK <- select(test_HK, Expression, Decay, Rate01Median, Rate10Median, EjectMedian, BurstMedian)


test_TF <- filter(Compare, Type == "TF"& rep == "rep1")

test_TF <- select(test_TF, Expression, Decay, Rate01Median, Rate10Median, EjectMedian, BurstMedian)

x <- test_HK[ , "Expression"]
y <- test_TF[ , "Expression"]


# Do x and y come from the same distribution?
ks.test(x, y)


##Test Expression rep2
test_HK <- filter(Compare, Type == "HK" & rep == "rep2")

test_HK <- select(test_HK, Expression, Decay, Rate01Median, Rate10Median, EjectMedian, BurstMedian)


test_TF <- filter(Compare, Type == "TF"& rep == "rep2")

test_TF <- select(test_TF, Expression, Decay, Rate01Median, Rate10Median, EjectMedian, BurstMedian)

x <- test_HK[ , "Expression"]
y <- test_TF[ , "Expression"]


# Do x and y come from the same distribution?
ks.test(x, y)



##Test Decay rep1
test_HK <- filter(Compare, Type == "HK" & rep == "rep1")

test_HK <- select(test_HK, Expression, Decay, Rate01Median, Rate10Median, EjectMedian, BurstMedian)


test_TF <- filter(Compare, Type == "TF"& rep == "rep1")

test_TF <- select(test_TF, Expression, Decay, Rate01Median, Rate10Median, EjectMedian, BurstMedian)

x <- test_HK[ , "Decay"]
y <- test_TF[ , "Decay"]


# Do x and y come from the same distribution?
ks.test(x, y)


##Test Decay rep2
test_HK <- filter(Compare, Type == "HK" & rep == "rep2")

test_HK <- select(test_HK, Expression, Decay, Rate01Median, Rate10Median, EjectMedian, BurstMedian)


test_TF <- filter(Compare, Type == "TF"& rep == "rep2")

test_TF <- select(test_TF, Expression, Decay, Rate01Median, Rate10Median, EjectMedian, BurstMedian)

x <- test_HK[ , "Decay"]
y <- test_TF[ , "Decay"]


# Do x and y come from the same distribution?
ks.test(x, y)



##Test Burst Size rep1
test_HK <- filter(Compare, Type == "HK" & rep == "rep1")

test_HK <- select(test_HK, Expression, Decay, Rate01Median, Rate10Median, EjectMedian, BurstMedian)


test_TF <- filter(Compare, Type == "TF"& rep == "rep1")

test_TF <- select(test_TF, Expression, Decay, Rate01Median, Rate10Median, EjectMedian, BurstMedian)

x <- test_HK[ , "BurstMedian"]
y <- test_TF[ , "BurstMedian"]


# Do x and y come from the same distribution?
ks.test(x, y)


##Test Burst Size rep2
test_HK <- filter(Compare, Type == "HK" & rep == "rep2")

test_HK <- select(test_HK, Expression, Decay, Rate01Median, Rate10Median, EjectMedian, BurstMedian)


test_TF <- filter(Compare, Type == "TF"& rep == "rep2")

test_TF <- select(test_TF, Expression, Decay, Rate01Median, Rate10Median, EjectMedian, BurstMedian)

x <- test_HK[ , "BurstMedian"]
y <- test_TF[ , "BurstMedian"]


# Do x and y come from the same distribution?
ks.test(x, y)



########

Compare <- transform(Compare, off_time = (1/Rate01Median))

###

##Test off_time rep1
test_HK <- filter(Compare, Type == "HK" & rep == "rep1")

test_HK <- select(test_HK, Expression, Decay, Rate01Median, Rate10Median, EjectMedian, BurstMedian, off_time)


test_TF <- filter(Compare, Type == "TF"& rep == "rep1")

test_TF <- select(test_TF, Expression, Decay, Rate01Median, Rate10Median, EjectMedian, BurstMedian, off_time)

x <- test_HK[ , "off_time"]
y <- test_TF[ , "off_time"]


# Do x and y come from the same distribution?
ks.test(x, y)


##Test off_time rep2
test_HK <- filter(Compare, Type == "HK" & rep == "rep2")

test_HK <- select(test_HK, Expression, Decay, Rate01Median, Rate10Median, EjectMedian, BurstMedian, off_time)


test_TF <- filter(Compare, Type == "TF"& rep == "rep2")

test_TF <- select(test_TF, Expression, Decay, Rate01Median, Rate10Median, EjectMedian, BurstMedian, off_time)

x <- test_HK[ , "off_time"]
y <- test_TF[ , "off_time"]


# Do x and y come from the same distribution?
ks.test(x, y)

##

Compare %>%
  group_by(Type, rep) %>%
  summarise_at(vars(Expression), list(name = mean))


Compare %>%
  group_by(Type, rep) %>%
  summarise_at(vars(off_time), list(name = median))

Compare %>%
  group_by(Type, rep) %>%
  summarise_at(vars(BurstMedian), list(name = median))


Compare%>%
  group_by(Type, rep) %>%
  summarise(
    offtime = median(BurstMedian),
    sd = sd(BurstMedian),
    n = n(),
    se = sd / sqrt(n)
  )

Compare%>%
  group_by(Type, rep) %>%
  summarise(
    offtime = median(off_time),
    sd = sd(off_time),
    n = n(),
    se = sd / sqrt(n)
  )



####Combine replicates - plots with error bars


ecdf_values_rep1_HK <- ecdf(Compare$Expression[Compare$rep == "rep1" & Compare$Type == "HK"])
ecdf_values_rep2_HK <- ecdf(Compare$Expression[Compare$rep == "rep2" & Compare$Type == "HK"])
ecdf_values_rep1_TF <- ecdf(Compare$Expression[Compare$rep == "rep1" & Compare$Type == "TF"])
ecdf_values_rep2_TF <- ecdf(Compare$Expression[Compare$rep == "rep2" & Compare$Type == "TF"])

Compare_rep1_HK <- data.frame(Expression = unique(Compare$Expression[Compare$Type == "HK"]), 
                              ECDF = ecdf_values_rep1_HK(unique(Compare$Expression[Compare$Type == "HK"])),
                              y = 1.96 * sqrt(ecdf_values_rep1_HK(unique(Compare$Expression[Compare$Type == "HK"]))) * 
                                sqrt(1 - ecdf_values_rep1_HK(unique(Compare$Expression[Compare$Type == "HK"]))) / 
                                sqrt(sum(Compare$rep == "rep1" & Compare$Type == "HK")))

Compare_rep1_HK$rep <- "rep1"
Compare_rep1_HK$Type <- "HK"

Compare_rep2_HK <- data.frame(Expression = unique(Compare$Expression[Compare$Type == "HK"]), 
                              ECDF = ecdf_values_rep2_HK(unique(Compare$Expression[Compare$Type == "HK"])),
                              y = 1.96 * sqrt(ecdf_values_rep2_HK(unique(Compare$Expression[Compare$Type == "HK"]))) * 
                                sqrt(1 - ecdf_values_rep2_HK(unique(Compare$Expression[Compare$Type == "HK"]))) / 
                                sqrt(sum(Compare$rep == "rep2" & Compare$Type == "HK")))
Compare_rep2_HK$rep <- "rep2"
Compare_rep2_HK$Type <- "HK"

Compare_rep1_TF <- data.frame(Expression = unique(Compare$Expression[Compare$Type == "TF"]), 
                              ECDF = ecdf_values_rep1_TF(unique(Compare$Expression[Compare$Type == "TF"])),
                              y = 1.96 * sqrt(ecdf_values_rep1_TF(unique(Compare$Expression[Compare$Type == "TF"]))) * 
                                sqrt(1 - ecdf_values_rep1_TF(unique(Compare$Expression[Compare$Type == "TF"]))) / 
                                sqrt(sum(Compare$rep == "rep1" & Compare$Type == "TF")))
Compare_rep1_TF$rep <- "rep1"
Compare_rep1_TF$Type <- "TF"

Compare_rep2_TF <- data.frame(Expression = unique(Compare$Expression[Compare$Type == "TF"]), 
                              ECDF = ecdf_values_rep2_TF(unique(Compare$Expression[Compare$Type == "TF"])),
                              y = 1.96 * sqrt(ecdf_values_rep2_TF(unique(Compare$Expression[Compare$Type == "TF"]))) * 
                                sqrt(1 - ecdf_values_rep2_TF(unique(Compare$Expression[Compare$Type == "TF"]))) / 
                                sqrt(sum(Compare$rep == "rep2" & Compare$Type == "TF")))
Compare_rep2_TF$rep <- "rep2"
Compare_rep2_TF$Type <- "TF"

combined_data <- rbind(Compare_rep1_HK, Compare_rep2_HK, Compare_rep1_TF, Compare_rep2_TF)

A <- ggplot(combined_data, aes(x = log10(Expression), y = ECDF, color = Type, linetype = rep)) +
  geom_step() +
  geom_errorbar(aes(ymin = ECDF - y, ymax = ECDF + y), width = 0) +
  xlab("log10 Expression [mean counts]") +
  ylab("fraction")  +
  theme(axis.text = element_text(size = 6), 
        axis.title = element_text(size = 7)) +   theme_classic2() + 
  scale_color_brewer(palette = "Set1")+
  guides(color = FALSE, linetype = FALSE)



#Decay
ecdf_values_rep1_HK <- ecdf(Compare$Decay[Compare$rep == "rep1" & Compare$Type == "HK"])
ecdf_values_rep2_HK <- ecdf(Compare$Decay[Compare$rep == "rep2" & Compare$Type == "HK"])
ecdf_values_rep1_TF <- ecdf(Compare$Decay[Compare$rep == "rep1" & Compare$Type == "TF"])
ecdf_values_rep2_TF <- ecdf(Compare$Decay[Compare$rep == "rep2" & Compare$Type == "TF"])

Compare_rep1_HK <- data.frame(Decay = unique(Compare$Decay[Compare$Type == "HK"]), 
                              ECDF = ecdf_values_rep1_HK(unique(Compare$Decay[Compare$Type == "HK"])),
                              y = 1.96 * sqrt(ecdf_values_rep1_HK(unique(Compare$Decay[Compare$Type == "HK"]))) * 
                               sqrt(1 - ecdf_values_rep1_HK(unique(Compare$Decay[Compare$Type == "HK"]))) / 
                               sqrt(sum(Compare$rep == "rep1" & Compare$Type == "HK")))


#######
Compare_rep1_HK$rep <- "rep1"
Compare_rep1_HK$Type <- "HK"

Compare_rep2_HK <- data.frame(Decay = unique(Compare$Decay[Compare$Type == "HK"]), 
                              ECDF = ecdf_values_rep2_HK(unique(Compare$Decay[Compare$Type == "HK"])),
y = 1.96 * sqrt(ecdf_values_rep2_HK(unique(Compare$Decay[Compare$Type == "HK"]))) * 
  sqrt(1 - ecdf_values_rep2_HK(unique(Compare$Decay[Compare$Type == "HK"]))) / 
  sqrt(sum(Compare$rep == "rep2" & Compare$Type == "HK")))
Compare_rep2_HK$rep <- "rep2"
Compare_rep2_HK$Type <- "HK"

Compare_rep1_TF <- data.frame(Decay = unique(Compare$Decay[Compare$Type == "TF"]), 
                              ECDF = ecdf_values_rep1_TF(unique(Compare$Decay[Compare$Type == "TF"])),
y = 1.96 * sqrt(ecdf_values_rep1_TF(unique(Compare$Decay[Compare$Type == "TF"]))) * 
  sqrt(1 - ecdf_values_rep1_TF(unique(Compare$Decay[Compare$Type == "TF"]))) / 
  sqrt(sum(Compare$rep == "rep1" & Compare$Type == "TF")))
Compare_rep1_TF$rep <- "rep1"
Compare_rep1_TF$Type <- "TF"

Compare_rep2_TF <- data.frame(Decay = unique(Compare$Decay[Compare$Type == "TF"]), 
                              ECDF = ecdf_values_rep2_TF(unique(Compare$Decay[Compare$Type == "TF"])),
y = 1.96 * sqrt(ecdf_values_rep2_TF(unique(Compare$Decay[Compare$Type == "TF"]))) * 
  sqrt(1 - ecdf_values_rep2_TF(unique(Compare$Decay[Compare$Type == "TF"]))) / 
  sqrt(sum(Compare$rep == "rep2" & Compare$Type == "TF")))
Compare_rep2_TF$rep <- "rep2"
Compare_rep2_TF$Type <- "TF"

combined_data <- rbind(Compare_rep1_HK, Compare_rep2_HK, Compare_rep1_TF, Compare_rep2_TF)

B <- ggplot(combined_data, aes(x = log10(Decay), y = ECDF, color = Type, linetype = rep)) +
  geom_step() +
  geom_errorbar(aes(ymin = ECDF - y, ymax = ECDF + y), width = 0) +
  xlab("log10 mRNA Decay [1/min]") +
  ylab("fraction")  +
  theme(axis.text = element_text(size = 6), 
        axis.title = element_text(size = 7)) +   theme_classic2() + 
  scale_color_brewer(palette = "Set1")+
  guides(color = FALSE, linetype = FALSE)




###off time


ecdf_values_rep1_HK <- ecdf(Compare$off_time[Compare$rep == "rep1" & Compare$Type == "HK"])
ecdf_values_rep2_HK <- ecdf(Compare$off_time[Compare$rep == "rep2" & Compare$Type == "HK"])
ecdf_values_rep1_TF <- ecdf(Compare$off_time[Compare$rep == "rep1" & Compare$Type == "TF"])
ecdf_values_rep2_TF <- ecdf(Compare$off_time[Compare$rep == "rep2" & Compare$Type == "TF"])

Compare_rep1_HK <- data.frame(off_time = unique(Compare$off_time[Compare$Type == "HK"]), 
                              ECDF = ecdf_values_rep1_HK(unique(Compare$off_time[Compare$Type == "HK"])),
                              y = 1.96 * sqrt(ecdf_values_rep1_HK(unique(Compare$off_time[Compare$Type == "HK"]))) * 
                                sqrt(1 - ecdf_values_rep1_HK(unique(Compare$off_time[Compare$Type == "HK"]))) / 
                                sqrt(sum(Compare$rep == "rep1" & Compare$Type == "HK")))


#######
Compare_rep1_HK$rep <- "rep1"
Compare_rep1_HK$Type <- "HK"

Compare_rep2_HK <- data.frame(off_time = unique(Compare$off_time[Compare$Type == "HK"]), 
                              ECDF = ecdf_values_rep2_HK(unique(Compare$off_time[Compare$Type == "HK"])),
                              y = 1.96 * sqrt(ecdf_values_rep2_HK(unique(Compare$off_time[Compare$Type == "HK"]))) * 
                                sqrt(1 - ecdf_values_rep2_HK(unique(Compare$off_time[Compare$Type == "HK"]))) / 
                                sqrt(sum(Compare$rep == "rep2" & Compare$Type == "HK")))
Compare_rep2_HK$rep <- "rep2"
Compare_rep2_HK$Type <- "HK"

Compare_rep1_TF <- data.frame(off_time = unique(Compare$off_time[Compare$Type == "TF"]), 
                              ECDF = ecdf_values_rep1_TF(unique(Compare$off_time[Compare$Type == "TF"])),
                              y = 1.96 * sqrt(ecdf_values_rep1_TF(unique(Compare$off_time[Compare$Type == "TF"]))) * 
                                sqrt(1 - ecdf_values_rep1_TF(unique(Compare$off_time[Compare$Type == "TF"]))) / 
                                sqrt(sum(Compare$rep == "rep1" & Compare$Type == "TF")))
Compare_rep1_TF$rep <- "rep1"
Compare_rep1_TF$Type <- "TF"

Compare_rep2_TF <- data.frame(off_time = unique(Compare$off_time[Compare$Type == "TF"]), 
                              ECDF = ecdf_values_rep2_TF(unique(Compare$off_time[Compare$Type == "TF"])),
                              y = 1.96 * sqrt(ecdf_values_rep2_TF(unique(Compare$off_time[Compare$Type == "TF"]))) * 
                                sqrt(1 - ecdf_values_rep2_TF(unique(Compare$off_time[Compare$Type == "TF"]))) / 
                                sqrt(sum(Compare$rep == "rep2" & Compare$Type == "TF")))
Compare_rep2_TF$rep <- "rep2"
Compare_rep2_TF$Type <- "TF"

combined_data <- rbind(Compare_rep1_HK, Compare_rep2_HK, Compare_rep1_TF, Compare_rep2_TF)

C <- ggplot(combined_data, aes(x = log10(off_time), y = ECDF, color = Type, linetype = rep)) +
  geom_step() +
  geom_errorbar(aes(ymin = ECDF - y, ymax = ECDF + y), width = 0) +
  xlab("log10 OFF time [min]") +
  ylab("fraction") +
  theme(axis.text = element_text(size = 6), 
        axis.title = element_text(size = 7)) +   theme_classic2() + 
  scale_color_brewer(palette = "Set1")+
  guides(color = FALSE, linetype = FALSE)


###########

#Burst Size

ecdf_values_rep1_HK <- ecdf(Compare$BurstMedian[Compare$rep == "rep1" & Compare$Type == "HK"])
ecdf_values_rep2_HK <- ecdf(Compare$BurstMedian[Compare$rep == "rep2" & Compare$Type == "HK"])
ecdf_values_rep1_TF <- ecdf(Compare$BurstMedian[Compare$rep == "rep1" & Compare$Type == "TF"])
ecdf_values_rep2_TF <- ecdf(Compare$BurstMedian[Compare$rep == "rep2" & Compare$Type == "TF"])

Compare_rep1_HK <- data.frame(BurstMedian = unique(Compare$BurstMedian[Compare$Type == "HK"]), 
                              ECDF = ecdf_values_rep1_HK(unique(Compare$BurstMedian[Compare$Type == "HK"])),
                              y = 2 * sqrt(ecdf_values_rep1_HK(unique(Compare$BurstMedian[Compare$Type == "HK"]))) * 
                                sqrt(1 - ecdf_values_rep1_HK(unique(Compare$BurstMedian[Compare$Type == "HK"]))) / 
                                sqrt(sum(Compare$rep == "rep1" & Compare$Type == "HK")))


#######
Compare_rep1_HK$rep <- "rep1"
Compare_rep1_HK$Type <- "HK"

Compare_rep2_HK <- data.frame(BurstMedian = unique(Compare$BurstMedian[Compare$Type == "HK"]), 
                              ECDF = ecdf_values_rep2_HK(unique(Compare$BurstMedian[Compare$Type == "HK"])),
                              y = 2 * sqrt(ecdf_values_rep2_HK(unique(Compare$BurstMedian[Compare$Type == "HK"]))) * 
                                sqrt(1 - ecdf_values_rep2_HK(unique(Compare$BurstMedian[Compare$Type == "HK"]))) / 
                                sqrt(sum(Compare$rep == "rep2" & Compare$Type == "HK")))
Compare_rep2_HK$rep <- "rep2"
Compare_rep2_HK$Type <- "HK"

Compare_rep1_TF <- data.frame(BurstMedian = unique(Compare$BurstMedian[Compare$Type == "TF"]), 
                              ECDF = ecdf_values_rep1_TF(unique(Compare$BurstMedian[Compare$Type == "TF"])),
                              y = 2 * sqrt(ecdf_values_rep1_TF(unique(Compare$BurstMedian[Compare$Type == "TF"]))) * 
                                sqrt(1 - ecdf_values_rep1_TF(unique(Compare$BurstMedian[Compare$Type == "TF"]))) / 
                                sqrt(sum(Compare$rep == "rep1" & Compare$Type == "TF")))
Compare_rep1_TF$rep <- "rep1"
Compare_rep1_TF$Type <- "TF"

Compare_rep2_TF <- data.frame(BurstMedian = unique(Compare$BurstMedian[Compare$Type == "TF"]), 
                              ECDF = ecdf_values_rep2_TF(unique(Compare$BurstMedian[Compare$Type == "TF"])),
                              y = 2 * sqrt(ecdf_values_rep2_TF(unique(Compare$BurstMedian[Compare$Type == "TF"]))) * 
                                sqrt(1 - ecdf_values_rep2_TF(unique(Compare$BurstMedian[Compare$Type == "TF"]))) / 
                                sqrt(sum(Compare$rep == "rep2" & Compare$Type == "TF")))
Compare_rep2_TF$rep <- "rep2"
Compare_rep2_TF$Type <- "TF"

combined_data <- rbind(Compare_rep1_HK, Compare_rep2_HK, Compare_rep1_TF, Compare_rep2_TF)

D <- ggplot(combined_data, aes(x = log10(BurstMedian), y = ECDF, color = Type, linetype = rep)) +
  geom_step() +
  geom_errorbar(aes(ymin = ECDF - y, ymax = ECDF + y), width = 0) +
  xlab("log10 Burst Size [mRNA/burst]") +
  ylab("fraction") +
  theme(axis.text = element_text(size = 6), 
        axis.title = element_text(size = 7)) +   theme_classic2() + 
  scale_color_brewer(palette = "Set1")+
  guides(color = FALSE, linetype = FALSE)

ggarrange(A, B, C, D, ncol = 2, nrow = 2)




