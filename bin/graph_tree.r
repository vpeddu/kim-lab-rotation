#library('ape')
library('ggtree')
library('viridis')
library('reshape2')
library('tidyverse')

args = commandArgs(trailingOnly=TRUE)



pancTree<-read.tree(args[1])
csvfilename = paste0(args[2],'_filtered.csv')
panc<-read.csv(csvfilename)


pancTree$tip.label<-sapply(strsplit(as.character(pancTree$tip.label), "_5__"), `[`, 1)
pancTree$tip.label<-sapply(strsplit(as.character(pancTree$tip.label), "k_"), `[`, 2)

#newlabels<-strsplit(pancTree$tip.label,'___')
#trim1<-unlist(newlabels)[2*(1:length(pancTree$tip.label)) - 1]
#trim2<-strsplit(trim1,'_5')
#trim3<-unlist(trim2)[2*(1:length(pancTree$tip.label)) - 1]
#pancTree$tip.label<-trim1

metadata<-panc[c('id', 'avg_count')]
metadata$id <- gsub(':', '__', metadata$id)
metadata<-metadata[metadata$id %in% pancTree$tip.label,]

#rownames(metadata)<-metadata$id
#metadata$id<-NULL


## ggtree patch from github (https://github.com/YuLab-SMU/ggtree/issues/315)
facet_widths <- function(p, widths) {
  if (!is.null(names(widths))) {
    ## if (is.ggtree(p) && !is.null(names(widths))) {
    ## .panel <- levels(p$data$.panel)
    .panel <- panel_col_levels(p)
    w <- rep(1, length=length(.panel))
    names(w) <- .panel
    w[names(widths)] <- widths
    widths <- w
  }
  gt <- ggplot_gtable(ggplot_build(p))
  for(i in seq_along(widths)) {
    j <- gt$layout$l[grep(paste0('panel-1-', i), gt$layout$name)]
    gt$widths[j] = widths[i] * gt$widths[j]
  }
  return(ggplotify::as.ggplot(gt))
}


filename = paste0(args[2], '_unique_alu_tree_plot.pdf')
metadata$subfamily<-substr(metadata$id, start = 4, stop = 4)
metadata$id<-metadata$id
f<-metadata
f<-melt(f)
print('here')
newtree<-ggtree(pancTree) %<+% metadata +    
  geom_tiplab(size=1.5, align=FALSE, linesize=.5) + 
  #ggtitle('panc Alus') + 
  theme_tree2() + 
  #geom_tiplab(offset = .6, hjust = .5) +
  geom_tippoint(aes(color = subfamily)) +
  theme(legend.position = "right") + scale_size_continuous(range = c(3, 10)) + 
  geom_facet(panel = "Average count", data = f, geom = ggstance::geom_barh,
            aes(x = avg_count, y = y), fill = "#39568CFF", 
           stat = "identity", limits = c(0,200)) +
  xlim_expand(c(0, 120), 'Average count') +
  scale_color_viridis(discrete=TRUE) 
 newtree + facet_widths(newtree, widths = c(3,1))

 #facet_widths(newtree, c('Average Count'  = .3))
pdf(file = filename,   width = 10,  height = 10) 
newtree
dev.off() 
