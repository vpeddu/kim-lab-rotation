# library('ape')
# library('ggtree')
# library('viridis')
# library('reshape2')
#library('tidyverse')

args = commandArgs(trailingOnly=TRUE)

covid_file<-read.csv(args[1])
covid<-covid_file[covid_file$log2FoldChange > 1 & covid_file$padj < 0.05 & covid_file$family == 'Alu',]
covid$id <- paste0(covid$gene,'_range=',covid$chr,':',covid$start,'-',covid$end)


panc_file<-read.csv(args[2])
panc<-panc_file[panc_file$log2FoldChange > 1 & panc_file$padj < 0.05 & panc_file$clade == 'SINE',]
panc$id <- paste0(panc$gene,'_range=',panc$chr,':',panc$start,'-',panc$end)


overlap <- panc$id[panc$id %in% covid$id]

pancUnique<-panc[! (panc$id %in% covid$id),]

covidUnique<-panc[! (covid$id %in% panc$id),]


write.table(covidUnique$id, file = "covidUnique.txt", sep = "\n",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)
write.table(pancUnique$id, file = "pancUnique.txt", sep = "\n",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

write.csv(file = 'covid_filtered.csv',covidUnique)
write.csv(file = 'panc_filtered.csv',pancUnique)
# cat(capture.output(print(covidUnique$id), file="covidUnique.txt"))
# cat(capture.output(print(pancUnique$id), file="pancUnique.txt"))

#### Run panc unique analysis
