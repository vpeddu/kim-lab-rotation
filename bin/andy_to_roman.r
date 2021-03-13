# extracting TEs from Andy's gut DE experiment 

gut_full <- read.csv('/Users/vikas/Documents/UCSC/rotations/Kim/te_de/andy_gut/pancreas.plasma.ev.long.RNA.normalized.deseq.gene.diseaseState_PDAC_vs_healthy.lfcShrink.csv')
roman_de <- read.csv('/Users/vikas/Documents/UCSC/rotations/Kim/te_de/panc_v_ctrl_te.ip.table.csv')
repeat_data <- read.csv('/Users/vikas/Documents/UCSC/rotations/Kim/te_de/repeat_bigboi.txt', sep = '\t')

#brekdown of repeat classes
table(repeat_data$repClass)

#extract just the sines
only_sine = repeat_data[repeat_data$repClass == 'SINE',]

#decrease search space
lf2c.1 <- gut_full[which( abs ( gut_full$lfcShrinkResults.log2FoldChange)> 1) ,] 

r <- c()
for (i in 1:nrow(lf2c.1)) {
  print(i)
  if( lf2c.1$X[i] %in% only_sine$repName){ 
    r<-append(r,i)
    }
  
}
