library('ape')
library('ggtree')
library('tidyverse')
library("ComplexHeatmap")

setwd('/Users/vikas/Documents/UCSC/rotations/Kim/kim-lab-rotation')

covid_file<-read.csv('/Users/vikas/Documents/UCSC/rotations/Kim/te_de/covid_v_ctrl_te.ip.table.csv')
covid<-covid_file[covid_file$log2FoldChange > 1 && covid_file$padj < 0.05 && covid_file$family == 'Alu',]
covid$uniqueName <- paste0(covid$gene,'_range=',covid$chr,':',covid$start,'-',covid$end)



panc_file<-read.csv('/Users/vikas/Documents/UCSC/rotations/Kim/te_de/panc_v_ctrl_te.ip.table.csv')
panc<-panc_file[panc_file$log2FoldChange > 1 && panc_file$padj < 0.05 && panc_file$family == 'Alu',]
panc$uniqueName <- paste0(panc$gene,'_range=',panc$chr,':',panc$start,'-',panc$end)


overlap <- panc$uniqueName[panc$uniqueName %in% covid$uniqueName]

covidTree<-read.tree('/Users/vikas/Documents/UCSC/rotations/Kim/kim-lab-rotation/covid_alu_subset alignment consensus tree.newick')
int<-as_tibble(covidTree)



newlabels<-strsplit(covidTree$tip.label,'k_')
trim1<-unlist(newlabels)[2*(1:length(covidTree$tip.label))]
trim2<-strsplit(trim1,'_5')
trim3<-unlist(trim2)[2*(1:length(covidTree$tip.label)) - 1]
covidTree$tip.label<-trim3

metadata<-covid[c('uniqueName', 'avg_count')]
metadata<-metadata[metadata$uniqueName %in% covidTree$tip.label,]

#rownames(metadata)<-metadata$uniqueName
#metadata$uniqueName<-NULL

newtree<-ggtree(covidTree) +    
  geom_tiplab(size=1.5, align=FALSE, linesize=.5) + 
  scale_x_ggtree() 
  #theme_tree2() 

# newtree+ 
#   geom_facet(panel = "Trait", data = metadata, geom=geom_segment)



p3 <- facet_plot(newtree, panel='Average Count', data=metadata, geom=geom_segment, 
                 aes(x=0, xend = avg_count, yend = y), size=3, color='blue4') 

p3 + facet_widths(p3, c('Average Count'  = .5))

ggsave(p3, height = 10, width = 10)

gheatmap(newtree, metadata, offset=5, width=0.5, font.size=3, 
         colnames_angle=-45, hjust=0) 
  # scale_fill_manual(breaks=c("HuH3N2", "pdm", "trig"), 
  #                   values=c("steelblue", "firebrick", "darkgreen"), name="genotype")
