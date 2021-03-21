setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(dplyr)
library(ggpubr)
library(RSvgDevice)

kmers <- read.csv("chm13.draft_v1.0.freq", header = FALSE, sep = "\t", na.strings="NA")

kmers <- kmers %>% mutate(missing = if_else(is.na(V2) | is.na(V3), 'missing', 'non_missing'))

kmers[is.na(kmers)] <- 0

g<-ggplot(kmers%>%arrange(desc(missing)), aes(x=V2, y=V3)) + 
  #stat_density_2d(aes(fill = ..level..), geom="polygon", n = 200) +
  geom_bin2d(bins=1000, size = 10) +
  scale_fill_gradient(low="#21409A", high="#ED1C24") +
  #  geom_smooth(method = lm, size = 5, color="black")+
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
    axis.ticks = element_line(colour = "black", size=2),
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
  xlab("Illumina K*")+
  ylab("Pacbio HiFi K*")+
  xlim(c(-25,25))+
  ylim(c(-25,25))+
  guides(fill = guide_colourbar(barwidth = 1.5, barheight = 100))+
  scale_size(range = c(0,20), breaks=c(100,1000,5000))
#scale_color_manual(breaks = c("missing", "non_missing"),
#                   values=c("red", "black"))

devSVG(file = "correlation_Kstar.svg", width = 50, height = 50, bg = "white", fg = "black",
       onefile = TRUE, xmlHeader = TRUE)
g
dev.off()

png("correlation_Kstar_4legend.png", width = 2000, height = 2000)
g
dev.off()
