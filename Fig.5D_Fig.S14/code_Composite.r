library(ggplot2)

annot <- data.frame(lab= c("-5kb", "Start", "End", "+5kb"), mpos = c(0,49,149,199))
lines <- data.frame( x = c(49, 149))

#JQ1 
data <- read.table("Proseq_JQ1_composite.txt", sep="\t", header=T, stringsAsFactors=F)
temp <- subset(data, gt=="DJQ1" | gt=="DDMSO")
ggplot(aes(pos,Mn, colour=gt), data=temp) +
  geom_line()+
  geom_vline(aes(xintercept = x), data = lines, linetype = 3) +
  scale_x_continuous(breaks = annot$mpos, labels = annot$lab)+
  ylim(c(0,7.5))+
  labs(x = "Position", y= "mean rpkm", title = "Proseq composite on each category for Down Sig genes")+
  theme_classic()+
  theme( axis.title = element_text(size = 7),  
         axis.text.y = element_text(size = 6))

temp <- subset(data, gt=="OJQ1" | gt=="ODMSO")
ggplot(aes(pos,Mn, colour=gt), data=temp) +
  geom_line()+
  geom_vline(aes(xintercept = x), data = lines, linetype = 3) +
  scale_x_continuous(breaks = annot$mpos, labels = annot$lab)+
  ylim(c(0,7))+
  labs(x = "Position", y= "mean rpkm", title = "Proseq composite on each category for Other genes")+
  theme_classic()+
  theme( axis.title = element_text(size = 7),  
         axis.text.y = element_text(size = 6))

##Supp
ggplot(aes(pos,Mn, colour=gt), data=temp) +
  geom_line()+
  geom_vline(aes(xintercept = x), data = lines, linetype = 3) +
  scale_x_continuous(breaks = annot$mpos, labels = annot$lab)+
  ylim(c(0.4,0.7))+
  labs(x = "Position", y= "mean rpkm", title = "Proseq composite on each category for Other genes")+
  theme_classic()+
  theme( axis.title = element_text(size = 7),  
         axis.text.y = element_text(size = 6))


#Med26Degron 
data <- read.table("Proseq_Med26Degron_composite.txt", sep="\t", header=T, stringsAsFactors=F)
temp <- subset(data, gt=="DIAA" | gt=="DDMSO")
ggplot(aes(pos,Mn, colour=gt), data=temp) +
  geom_line()+
  geom_vline(aes(xintercept = x), data = lines, linetype = 3) +
  scale_x_continuous(breaks = annot$mpos, labels = annot$lab)+
  ylim(c(0,7.5))+
  labs(x = "Position", y= "mean rpkm", title = "Proseq composite on each category for Down Sig genes")+
  theme_classic()+
  theme( axis.title = element_text(size = 7),  
         axis.text.y = element_text(size = 6))

temp <- subset(data, gt=="OIAA" | gt=="ODMSO")
ggplot(aes(pos,Mn, colour=gt), data=temp) +
  geom_line()+
  geom_vline(aes(xintercept = x), data = lines, linetype = 3) +
  scale_x_continuous(breaks = annot$mpos, labels = annot$lab)+
  ylim(c(0,7))+
  labs(x = "Position", y= "mean rpkm", title = "Proseq composite on each category for Other genes")+
  theme_classic()+
  theme( axis.title = element_text(size = 7),  
         axis.text.y = element_text(size = 6))

## Supp

ggplot(aes(pos,Mn, colour=gt), data=temp) +
  geom_line()+
  geom_vline(aes(xintercept = x), data = lines, linetype = 3) +
  scale_x_continuous(breaks = annot$mpos, labels = annot$lab)+
  ylim(c(0.4,0.7))+
  labs(x = "Position", y= "mean rpkm", title = "Proseq composite on each category for Other genes")+
  theme_classic()+
  theme( axis.title = element_text(size = 7),  
         axis.text.y = element_text(size = 6))

