library('ape')
library('ggtree')
library('viridis')
library('reshape2')
library('tidyverse')
library("ComplexHeatmap")



setwd('/Users/vikas/Documents/UCSC/rotations/Kim/kim-lab-rotation/unique_plot')

covid_file<-read.csv('/Users/vikas/Documents/UCSC/rotations/Kim/te_de/covid_v_ctrl_te.ip.table.csv')
covid<-covid_file[covid_file$log2FoldChange > 1 & covid_file$padj < 0.05 & covid_file$family == 'Alu',]
covid$id <- paste0(covid$gene,'_range=',covid$chr,':',covid$start,'-',covid$end)


panc_file<-read.csv('/Users/vikas/Documents/UCSC/rotations/Kim/te_de/panc_v_ctrl_te.ip.table.csv')
panc<-panc_file[panc_file$log2FoldChange > 1 & panc_file$padj < 0.05 & panc_file$clade == 'SINE',]
panc$id <- paste0(panc$gene,'_range=',panc$chr,':',panc$start,'-',panc$end)


overlap <- panc$id[panc$id %in% covid$id]

pancUnique<-panc[! (panc$id %in% covid$id),]

covidUnique<-panc[! (covid$id %in% panc$id),]


#### Run panc unique analysis

# grab fasta records
print('getting fastas')
for (i in 1:nrow(pancUnique)){ 
  grab = paste0('grep -A1 "', pancUnique$id[i], '" unwrapped.fasta >> panc_unique_alu_subset.fasta ')
  system(grab)
}
# mafft align
print('aligning')
  align = paste0('MAFFT panc_unique_alu_subset.fasta > panc_unique_alu_aligned.fasta')
  system(align)
# fasttree? 
print('building tree')
  buildtree = paste0('FastTree -gtr -nt < panc_unique_alu_aligned.fasta > panc_unique_alu_tree.newick')


pancTree<-read.tree('panc_unique_alu_tree.newick')

newlabels<-strsplit(pancTree$tip.label,'k_')
trim1<-unlist(newlabels)[2*(1:length(pancTree$tip.label))]
trim2<-strsplit(trim1,'_5')
trim3<-unlist(trim2)[2*(1:length(pancTree$tip.label)) - 1]
pancTree$tip.label<-trim3

metadata<-panc[c('id', 'avg_count')]
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



metadata$subfamily<-substr(metadata$id, start = 4, stop = 4)
metadata$id<-metadata$id
f<-metadata
f<-melt(f)
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
 newtree
 pdf(file = "panc_unique_alu_tree_plot.pdf",   width = 4,  height = 4) 
# 
# f = metadata
#  facet_plot(newtree,inherit.aes= TRUE,  panel='Average Count',
#                   geom=geom_segment,
#                   data = metadata,
#                   aes(x= 0, xend = avg_count, y=y, yend=y),
#                   color = '#39568CFF', size=3) 
#  p3
# facet_widths(p3, c('Average Count'  = .3))






# 
# gheatmap(newtree, metadata, offset=5, width=0.5, font.size=3, 
#          colnames_angle=-45, hjust=0) 
  # scale_fill_manual(breaks=c("HuH3N2", "pdm", "trig"), 
  #                   values=c("steelblue", "firebrick", "darkgreen"), name="genotype")
