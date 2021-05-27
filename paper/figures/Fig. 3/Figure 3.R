setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)

text_size1=7
text_size2=80

#################################
#QV and QV*
#################################

AppendMe <- function(dfNames) {
  do.call(rbind, lapply(dfNames, function(x) {
    cbind(get(x), source = x)
  }))
}

species_fun <- function(variable,value){
  return(species_titles[value])
}

grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(
    legend.position="bottom", 
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.key = element_rect(fill = NA, color = NA)))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(arrangeGrob(
    arrangeGrob(plots[[1]] + theme(legend.position="none", plot.margin = unit(c(0,2,2,8), "cm")),
                top = textGrob("a", x = unit(0, "npc"), y = unit(0.3, "npc"), just=c("left","top"),
                               gp=gpar(fontsize=rel(text_size2), fontface="bold"))),
    arrangeGrob(plots[[2]] + theme(legend.position="none", plot.margin = unit(c(0,2,2,0), "cm")),
                top = textGrob("b", x = unit(-0.1, "npc"), y = unit(0.3, "npc"), just=c("left","top"),
                               gp=gpar(fontsize=rel(text_size2), fontface="bold"))),
    arrangeGrob(plots[[3]] + theme(legend.position="none", plot.margin = unit(c(0,0,2,0), "cm")),
                top = textGrob("c", x = unit(-0.1, "npc"), y = unit(0.3, "npc"), just=c("left","top"),
                               gp=gpar(fontsize=rel(text_size2), fontface="bold"))),
    nrow=1),
               legend, nrow=2, heights = unit.c(unit(1, "npc") - lheight, lheight))
}

species_titles<-list('fArcCen1' = "fArcCen1 (VGP v1.0)",
                     'rGopEvg1' = "rGopEvg1 (VGP v1.5)",
                     'bTaeGut1' = "bTaeGut1 (VGP v1.6)")

plot_lines <- function(df, stage, qv, dataset, title, y_title, steps) {
  df<-df %>% filter(type==dataset) %>% filter(stage %in% steps)
  plt <- ggplot(df, aes(x=stage, y=qv)) +
    geom_point(size=8, aes(color=as.factor(legend)))+
    geom_line(size=2, aes(group=as.factor(merfin), color=as.factor(legend), linetype=as.factor(legend)))+
    labs(x="Stage", y = y_title)+
    theme(
      panel.background = element_rect(fill = "white",
                                      colour = "gray",
                                      size = 2, linetype = "solid"),
      panel.grid.major = element_line(size = 1, linetype = 'dotted',
                                      colour = "grey"), 
      panel.grid.minor = element_line(size = 1, linetype = 'dotted',
                                      colour = "grey"),
      plot.title = element_text(family = "Arial", hjust = 0.5,
                                size = rel(text_size1), colour = "black", face="bold",
                                margin = margin(t = 10, r = 0, b = 10, l = 0)), 
      axis.title = element_text(family = "Arial",
                                size = rel(text_size1), colour = "black", face="bold"), 
      axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 40)),
      axis.title.x = element_blank(),
      axis.text = element_text(family = "Arial",
                               size = rel(text_size1), colour = "black"),
      axis.ticks = element_line(size = 2),
      axis.ticks.length = unit(0.5,"cm"),
      axis.text.x=element_text(margin = margin(t = 20), angle = 45, hjust = 1),
      axis.text.y.right=element_text(margin = margin(l = 10)),
      plot.tag = element_text(size = rel(text_size1), face = "bold"),
      legend.text = element_text(size = rel(text_size1), face = "bold"),
      legend.position="bottom",
      legend.key.size = unit(2,"cm"),
      legend.title = element_blank(),
      legend.background = element_blank(),
      legend.key = element_rect(fill = NA, color = NA),
      legend.spacing.x = unit(1.0, 'cm'),
      strip.background = element_blank(), strip.text = element_blank(),
      panel.spacing = unit(2, "lines")
      )+
    ggtitle(title)+
    scale_color_manual(values=c("black","#3DDC97","#FF495C"), breaks = c(2,3,4), labels = rev(c("Default", "Merfin", "Merfin + Default")), guide = guide_legend(reverse = TRUE))+
    scale_linetype_manual(values=c("dashed","solid","solid"), breaks = c(2,3,4), labels = rev(c("Default", "Merfin", "Merfin + Default")), guide = guide_legend(reverse = TRUE))+
    scale_x_discrete(labels = c('Unpolished','Arrow','Freebayes 1','Freebayes 2'), expand = c(0.05,0.05))+
    scale_y_continuous(position = "right")
  plt + facet_grid(species ~ ., scales = "free", switch="both",
                   labeller = species_fun)
}

data <- read.csv("data.qv", header = TRUE, sep = "\t", na.strings="NA")
data$species <- factor(data$species, levels = c("fArcCen1", "rGopEvg1", "bTaeGut1"))

gr1<-plot_lines(data, 'stage', 'qv', 'merqury',
                     'Primary/Alternate combined',
                     'Merqury QV (missing kmers only)',
                     c('s4','t1','t2','t3'))

gr2<-plot_lines(data, 'stage', 'qv', 'merqury_lookup',
                     'Primary/Alternate combined',
                     'Adjusted Merqury QV (low frequency kmers)',
                     c('s4','t1','t2','t3'))

gr3<-plot_lines(data, 'stage', 'qv', 'lookup',
                     'Primary/Alternate combined',
                     'QV* (overrepresented kmers included)',
                     c('s4','t1','t2','t3'))

qv<-grid_arrange_shared_legend(gr1,gr2,gr3)

#################################
#gene annotation
#################################

grid_arrange_shared_legend2 <- function(plots) {
  g <- ggplotGrob(plots[[1]] + theme(
    legend.position="bottom", 
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.key = element_rect(fill = NA, color = NA)))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(arrangeGrob(
    arrangeGrob(plots[[1]] + theme(legend.position="none", plot.margin = unit(c(0,2,2,0), "cm")),
                top = textGrob("d", x = unit(0, "npc"), y = unit(0.3, "npc"), just=c("left","top"),
                               gp=gpar(fontsize=rel(text_size2), fontface="bold"))),
    arrangeGrob(plots[[2]] + theme(legend.position="none", plot.margin = unit(c(0,2,2,0), "cm")),
                top = textGrob("e", x = unit(0, "npc"), y = unit(0.3, "npc"), just=c("left","top"),
                               gp=gpar(fontsize=rel(text_size2), fontface="bold"))),
    arrangeGrob(plots[[3]] + theme(legend.position="none", plot.margin = unit(c(0,2,2,0), "cm")),
                top = textGrob("f", x = unit(0, "npc"), y = unit(0.3, "npc"), just=c("left","top"),
                               gp=gpar(fontsize=rel(text_size2), fontface="bold"))),
    nrow=1),
    legend, nrow=2, heights = unit.c(unit(1, "npc") - lheight, lheight))
}

plot_bars <- function(df, stage, readout, title) {
  plt<-ggplot(df, aes_string(x="stage", y=readout)) +
    geom_bar(stat="identity", aes(fill=stage), color="black",size=3)+
    theme(
    panel.background = element_rect(fill = "white",
                                    colour = "gray",
                                    size = 2, linetype = "solid"),
    panel.grid.major = element_line(size = 1, linetype = 'dotted',
                                    colour = "grey"), 
    panel.grid.minor = element_line(size = 1, linetype = 'dotted',
                                    colour = "grey"),
    plot.title = element_text(family = "Arial", hjust = 0.5,
                              size = rel(text_size1), colour = "black", face="bold",
                              margin = margin(t = 10, r = 0, b = 10, l = 0)), 
    axis.title = element_blank(),
    axis.text = element_text(family = "Arial",
                             size = rel(text_size1), colour = "black"),
    axis.ticks = element_line(size = 2),
    axis.ticks.length = unit(0.5,"cm"),
    axis.text.x=element_text(margin = margin(t = 20), angle = 45, hjust = 1),
    axis.text.y.right=element_text(margin = margin(l = 10)),
    plot.tag = element_text(size = rel(text_size1), face = "bold"),
    legend.text = element_text(size = rel(text_size1), face = "bold"),
    legend.position="bottom",
    legend.key.size = unit(2,"cm"),
    legend.title = element_blank(),
    legend.background = element_blank(),
    legend.key = element_rect(fill = NA, color = NA),
    legend.spacing.x = unit(1.0, 'cm'),
    strip.background = element_blank(), strip.text = element_blank(),
    panel.spacing = unit(2, "lines")
  )+
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
    scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
    ggtitle(title) +
    facet_grid(species ~ ., scales = "free", switch="both")
}

data2 <- read.csv("data.annotation", header = TRUE, sep = "\t", na.strings="NA")
data2$species <- factor(data2$species, levels = c("fArcCen1", "rGopEvg1"))
data2$stage <- factor(data2$stage, levels = c("Unpolished", "Polished", "Merfin"))

gr1<-plot_bars(data2, 'stage', 'codebreaks', 'Premature stop')

gr2<-plot_bars(data2, 'stage', 'frameshifts', 'Frameshifts')

gr3<-plot_bars(data2, 'stage', 'low.quality', 'Low quality')

#####Plotting

png("Figure 3 top.png", width = 6500, height = 5000)

annotation<-grid_arrange_shared_legend2(list(gr1,gr2,gr3))

grid.arrange(qv,annotation, nrow=1, widths=c(2.3,1))

dev.off()

