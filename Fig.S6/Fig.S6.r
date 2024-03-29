library("ggplot2")

dataP<- read.table("Brd4_JQ1.txt", sep="\t", header=T, stringsAsFactors=F)
dataP$FC <- as.numeric(dataP$FC)

ggplot(dataP,aes(cat,log2(FC)))+
  geom_violin()+
  geom_boxplot(width=.1, position=position_dodge(1))+
  geom_hline(yintercept=0, lty=3, color="black")+
  geom_hline(yintercept=median(log2(dataP[dataP$cat=="promoter",]$FC)), lty=3, color="blue")+
  geom_hline(yintercept=median(log2(dataP[dataP$cat=="SE",]$FC)), lty=3, color="red")+
  labs(x='log2(FC (JQ/DMSO))', title="Brd4 signal intensity FC ")

