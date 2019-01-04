rm(list=ls())
.libPaths(c('/Users/ssapkot/Documents/SAP/GenSel/MyRlibs',.libPaths()))

library(pegas)
setwd("/Users/ssapkot/Documents/GenomicSelection1/SAP_Yield/Additional_analysis/")
hapmap <- read.table("/Users/ssapkot/Documents/GenomicSelection1/SAP_Yield/01_Data/SAP_GS_caudatum_all.hmp.txt",row.names=1, sep="\t", stringsAsFactors = FALSE)
hmp <- as.data.frame(t(hapmap))
hmp <- hmp[-c(1:10),]
sites <- as.data.frame(colnames(hmp))
####pegas###heterozygosity for all sites
Result <- c()
for (i in 2:ncol(hmp)) {
  het <- c()
  h <- hmp[,i]
  snp <- as.character(sites[i,])
  het <- heterozygosity(h,variance=TRUE)
  het <- as.numeric(het)
  Result <- rbind(Result,c(snp,het))
}
colnames(Result) <- c("Site","Het","Var")
Result_Total <- as.data.frame(Result)
write.csv(Result_Total,"Heterzygosity_bysites_Caudatum.csv")
#####

