library(ggplot2)

data <- read.table("r220913_Med26Degron_log10copies.txt", sep="\t", header=T, stringsAsFactors=F)

data$c <- densCols(data$NTC, data$sgRNA, colramp = colorRampPalette(c("steelblue4", "grey80")))
data$c[data$is.spike] <- "red2"
ggplot(data) +
    geom_point(aes(NTC, sgRNA, color = c, size = is.spike)) +
    geom_abline(intercept = 0, slope = 1, col = "grey40") +
    scale_color_identity() +
    scale_size_manual("", values = c("TRUE" = 3, "FALSE" = 2)) +
    scale_x_continuous("DMSO [copies / cell]", breaks = -1:3, labels = 10^(-1:3),
                       limits = c(-1, 3)) +
    scale_y_continuous("-MED26 [copies / cell]", breaks = -1:3, labels = 10^(-1:3), limits = c(-1, 3)) +

    coord_equal(ratio = 1) +
    theme_bw() +
    theme(legend.position = "none", 
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())

 
