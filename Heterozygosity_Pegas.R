
library(pegas)

## read and format the hapmap file
hapmap <- read.table("01_Data/SAP_GS_caudatum_all.hmp.txt",row.names=1, sep="\t", stringsAsFactors = FALSE)
hmp <- as.data.frame(t(hapmap))
hmp <- hmp[-c(1:10),]

## extract the names of the snp sites
sites <- as.data.frame(colnames(hmp))

## heterozygosity for all sites
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

