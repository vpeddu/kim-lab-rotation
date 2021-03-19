library(ggplot2)
library(ggpol)
library(viridis)

setwd('/Users/vikas/Documents/UCSC/rotations/Kim/kim-lab-rotation')

# read panc data 
load('/Users/vikas/Documents/UCSC/rotations/Kim/kim-lab-rotation/alu_output/treeFigures/panc2021-02-19 23:48:02.rds')
pancgenes<-metadata
pancgenes$type<-'Pancreatic adenocarcinoma'

# read covid alu data 
load('/Users/vikas/Documents/UCSC/rotations/Kim/kim-lab-rotation/alu_output/treeFigures/covid2021-02-19 23:47:40.rds')
covidalugenes<-metadata
covidalugenes$type<-'Covid-19'

# read covid line data
load('/Users/vikas/Documents/UCSC/rotations/Kim/kim-lab-rotation/line_output/treeFigures/covid2021-02-22 21:24:11.rds')
covidlinegenes<-metadata
covidlinegenes$type<-'Covid-19'

alus<-rbind(pancgenes, covidalugenes)
line<-covidlinegenes
alu_boxplot<-ggplot(alus, aes(x = type, y = avg_count, color = subfamily)) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.4) + 
  geom_jitter(shape=16, position=position_jitter(0.2), size = 2) +
  scale_y_log10() +
  ylim(1,250) + 
  theme_classic() + 
  theme(legend.title = element_blank(), axis.title.x = element_blank()) + 
  ylab('CPM') + 
  scale_color_viridis_d()
alu_boxplot
ggsave(plot =alu_boxplot, 'alu_boxplot.pdf', height = 5, width = 5)

line_boxplot<-ggplot(line, aes(x = type, y = avg_count, color = subfamily)) + 
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=.4) + 
  geom_jitter(shape=16, position=position_jitter(0.2), size = 2) +
  scale_y_log10() +
  ylim(1,60) + 
  theme_classic() + 
  theme(legend.title = element_blank(), axis.title.x = element_blank()) + 
  ylab('CPM') + 
  scale_color_viridis_d()
line_boxplot
ggsave(plot =line_boxplot, 'line_boxplot.pdf', height = 5, width = 5)
