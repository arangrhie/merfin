setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(dplyr)
library(ggpubr)

kmers <- read.csv("illumina_vs_hifi20kb.txt", header = FALSE, sep = "\t", na.strings="NA")

kmers <- kmers %>% mutate(missing = if_else(is.na(V2) | is.na(V3), 'missing', 'non_missing'))

png("correlation_readK.png", width = 2000, height = 2000)

# Add regression line
ggplot(kmers, aes(x=V2, y=V3)) + geom_point() + geom_smooth(method = lm)+
  theme(
    plot.margin = margin(t = 20, b = 20),
    axis.title.x = element_text(family = "Arial",
                                size = rel(4),colour = "black"), 
    axis.title.y = element_text(family = "Arial",
                                size = rel(4),colour = "black"), 
    axis.text = element_text(family = "Arial",
                             size = rel(4), colour = "black"),
    axis.ticks.length = unit(.50, "cm"),
    axis.ticks = element_line(colour = "black", size=(0.6)),
    legend.text=element_text(size=30),
    legend.title=element_blank()
  )+
  xlab("Illumina PCR free readK")+
  ylab("Pacbio Hifi 20kb K*")+
  stat_cor(method = "pearson")+
  xlim(c(0,250000))+
  ylim(c(0,250000))


dev.off()

kmers <- read.csv("illumina_vs_hifi20kb.txt", header = FALSE, sep = "\t", na.strings="NA")

kmers <- kmers %>% mutate(missing = if_else(is.na(V2) | is.na(V3), 'missing', 'non_missing'))

kmers[is.na(kmers)] <- 0

png("correlation_Kstar.png", width = 5000, height = 5000)

# Add regression line
b<-ggplot(kmers%>%arrange(desc(missing)), aes(x=V2, y=V3)) + 
  #stat_density_2d(aes(fill = ..level..), geom="polygon", n = 200) +
  geom_bin2d(bins=1000, size = 10) +
  scale_fill_gradient(low="blue", high="red") +
  geom_smooth(method = lm, size = 5, color="black")+
  geom_point(data = kmers%>%filter(missing == "missing"), aes(size = V1), shape=5, stroke = 2) +
  theme(
    plot.margin = margin(t = 20, b = 20, l = 20),
    axis.title.x = element_text(family = "Arial",
                                size = rel(8),colour = "black"), 
    axis.title.y = element_text(family = "Arial",
                                size = rel(8),colour = "black"), 
    axis.text = element_text(family = "Arial",
                             size = rel(8), colour = "black"),
    axis.ticks.length = unit(1, "cm"),
    axis.ticks = element_line(colour = "black", size=(2)),
    legend.text=element_text(size=100),
    legend.title=element_blank(),
    legend.key = element_rect(fill = NA, color = NA),
    panel.background = element_rect(fill = "white",
                                    colour = "white",
                                    size = 2, linetype = "solid"),
    panel.grid.major = element_line(size = 1, linetype = 'solid',
                                    colour = "grey"), 
    panel.grid.minor = element_line(size = 1, linetype = 'solid',
                                    colour = "grey")
  )+
  xlab("Illumina PCR free K*")+
  ylab("Pacbio Hifi 20kb K*")+
  xlim(c(-25,25))+
  ylim(c(-25,25))+
  guides(fill = guide_colourbar(barwidth = 1.5, barheight = 100))+
  scale_size(range = c(5,20))
  #scale_color_manual(breaks = c("missing", "non_missing"),
  #                   values=c("red", "black"))

b

dev.off()

# Point + regression line
# Remove the confidence interval 
b + geom_point() + 
  geom_smooth(method = lm, se = FALSE)

# loess method: local regression fitting
b + geom_point() + geom_smooth()