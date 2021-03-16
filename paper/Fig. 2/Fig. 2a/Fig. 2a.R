setwd(dirname(rstudioapi::getSourceEditorContext()$path))

library(ggplot2)

#################################
#Fig 2a - Density plot
#################################

AppendMe <- function(dfNames) {
  do.call(rbind, lapply(dfNames, function(x) {
    cbind(get(x), source = x)
  }))
}

plot_multi_histogram <- function(df, feature, label_column, adjust, xlim, scale, labels) {
  plt <- ggplot(df, aes(x=eval(parse(text=feature)), color=eval(parse(text=label_column)))) +
    geom_density(size=4, adjust = adjust, aes(weight=frequency), show.legend=FALSE)+
    stat_density(aes(weight=frequency),
                 geom="line",position="identity", size = 0)+
    guides(colour = guide_legend(override.aes=list(size=4)))+
    labs(x="K*", y = "Density")+
    theme(
      panel.background = element_rect(fill = "white",
                                      colour = "white",
                                      size = 1, linetype = "solid"),
      panel.grid.major = element_line(size = 1, linetype = 'solid',
                                      colour = "grey"), 
      panel.grid.minor = element_line(size = 1, linetype = 'solid',
                                      colour = "grey"),
      axis.title = element_text(family = "Arial",
                                size = rel(5), colour = "black"), 
      axis.text = element_text(family = "Arial",
                               size = rel(5), colour = "black"),
      axis.ticks.length = unit(2, "cm"),
      axis.ticks = element_line(colour = "black", size=(1)),
      plot.tag = element_text(size = rel(5), face = "bold"),
      legend.text = element_text(size = rel(5)),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key = element_rect(fill = NA, color = NA),
      legend.key.width = unit(3,"cm")
    )
  if (!missing(labels)){
    plt <- plt + scale_color_discrete(labels = labels)
  }
  if (!missing(scale)){
    plt <- plt + scale_y_continuous(trans=scale)
  }
  if (missing(xlim)){
    plt
  }else{
    plt +
      xlim(xlim)
  }
}

t2t_chm13.20200602.hifi <- read.csv("hist/t2t-chm13.20200602.hifi", header = FALSE, sep = "\t", col.names = c("k", "frequency"), na.strings="NA")
t2t_chm13.20200611.hifi <- read.csv("hist/t2t-chm13.20200611.hifi", header = FALSE, sep = "\t", col.names = c("k", "frequency"), na.strings="NA")
t2t_chm13.20200727.hifi <- read.csv("hist/t2t-chm13.20200727.hifi", header = FALSE, sep = "\t", col.names = c("k", "frequency"), na.strings="NA")
t2t_chm13.20200904.hifi <- read.csv("hist/t2t-chm13.20200904.hifi", header = FALSE, sep = "\t", col.names = c("k", "frequency"), na.strings="NA")
t2t_chm13.v1.0.hifi <- read.csv("hist/chm13.draft_v1.0.hifi", header = FALSE, sep = "\t", col.names = c("k", "frequency"), na.strings="NA")
GCA_009914755.1_T2T_v0.7_genomic.hifi <- read.csv("hist/GCA_009914755.1_T2T_v0.7_genomic.hifi", header = FALSE, sep = "\t", col.names = c("k", "frequency"), na.strings="NA")


t2t_chm13.20200602.illumina <- read.csv("hist/t2t-chm13.20200602.illumina", header = FALSE, sep = "\t", col.names = c("k", "frequency"), na.strings="NA")
t2t_chm13.20200611.illumina <- read.csv("hist/t2t-chm13.20200611.illumina", header = FALSE, sep = "\t", col.names = c("k", "frequency"), na.strings="NA")
t2t_chm13.20200727.illumina <- read.csv("hist/t2t-chm13.20200727.illumina", header = FALSE, sep = "\t", col.names = c("k", "frequency"), na.strings="NA")
t2t_chm13.20200904.illumina <- read.csv("hist/t2t-chm13.20200904.illumina", header = FALSE, sep = "\t", col.names = c("k", "frequency"), na.strings="NA")
t2t_chm13.20200918.illumina <- read.csv("hist/t2t-chm13.20200918.illumina", header = FALSE, sep = "\t", col.names = c("k", "frequency"), na.strings="NA")
t2t_chm13.v1.0.illumina <- read.csv("hist/chm13.draft_v1.0.illumina", header = FALSE, sep = "\t", col.names = c("k", "frequency"), na.strings="NA")
GCA_009914755.1_T2T_v0.7_genomic.illumina <- read.csv("hist/GCA_009914755.1_T2T_v0.7_genomic.illumina", header = FALSE, sep = "\t", col.names = c("k", "frequency"), na.strings="NA")

hifi_df<-AppendMe(c("GCA_009914755.1_T2T_v0.7_genomic.hifi","t2t_chm13.v1.0.hifi"))

illumina_df<-AppendMe(c("GCA_009914755.1_T2T_v0.7_genomic.illumina","t2t_chm13.v1.0.illumina"))

png("hist_hifi.png", width = 2000, height = 2000)
plot_multi_histogram(hifi_df, 'k', 'source', 1, c(-150,150), labels = c('Hicanu assembly', 'T2T v1.0'))
dev.off()

png("hist_illumina.png", width = 2000, height = 2000)
plot_multi_histogram(illumina_df, 'k', 'source', 1, c(-150,150), labels = c('Hicanu assembly', 'T2T v1.0'))
dev.off()
