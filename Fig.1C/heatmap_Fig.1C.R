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
heatmap <- read.csv(file = 'FIISH_scRNA_heatmap.csv')



ggplot(heatmap, aes(x = Rate, y = Expression, fill = Value)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "red", limits = c(0.5, 1)) +
  coord_fixed() +theme_classic2()  + 
  geom_text(aes(Rate, Expression, label = Value), color = "white", size = 6) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(-2, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))


