library(readxl)

setwd('/Users/vikas/Documents/UCSC/rotations/Kim/kim-lab-rotation')
covid <- readxl::read_excel('/Users/vikas/Documents/UCSC/rotations/Kim/te_de/significant_subset.xlsx', sheet = 1)

covid<-covid[covid$family == 'Alu',]

covid$newname <- paste0(covid$gene,'_range=',covid$chr,':',covid$start,'-',covid$end)

for (i in 1:nrow(covid)){ 
  cmd = paste0('grep -A1 "', covid$newname[i], '" unwrapped.fasta >> covid_alu_subset.fasta ')
  system(cmd)
}
