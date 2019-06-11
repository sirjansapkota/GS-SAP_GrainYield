
.libPaths(c('/anaconda2/lib/R/library',.libPaths()))

library(rrBLUP)
library(ggplot2)

my.read.vcf <- function(file, special.char="##", ...) {
  my.search.term <- paste0(special.char, ".*")  # Making a search term that looks like: "##.*", tells R to find anything containing the pattern "##" followed by anything (* is wildcard)
  clean.lines <- sub(my.search.term, "", readLines(file)) # Replace any line containing the search term with nothing (in other words remove it)
  clean.lines2 <- sub("#CHROM", "CHROM", clean.lines) # Replace the #CHROM term in the header with CHROM, so R doesn't treat it as a special character
  read.table(..., text=paste(clean.lines2, collapse="\n")) # Pass the cleaned up lines to read.table
}

setwd("/Users/ssapkot/Documents/Experiments/SAP_GS_PopStr")
Y <- read.csv('data/BLUEs_pheno_all.csv', header=T, row.names = 1)
Y$SN = 1:nrow(Y)


##Setting up a loop for running for each cluster
for (a in 1:5){
    
    Y1 = Y[Y$Cluster == a,]
    nr = floor(nrow(Y1)/5)
    cvf <- data.frame(matrix("",nrow=nr*5, ncol=0))
    
    for (b in 1:100){
        
        c.list = vector("list", 5)
        
        ctemp = sample(Y1$SN,5*nr)
        c.list[[b]] = ctemp
        cvf[b] <- cbind(c.list[[b]])
    }
    fold <- rep(1:5, times = nr)
    cvf <- cbind(fold,cvf)

    print(nrow(cvf))
    }


GBS=my.read.vcf(file= "data/SAP_all_taxa.vcf", header=TRUE, stringsAsFactors = TRUE, as.is=TRUE)


f.column <- grep("FORMAT", colnames(GBS))

##Function Parse vcf file to convert to -1,0,1 format
parse.GBS <- function(x) {
  unique.x <- unique(x)
  alleles <- setdiff(unique.x,union("H","N"))
  y <- rep(0,length(x))
  y[which(x==alleles[1])] <- -1
  y[which(x==alleles[2])] <- 1
  y[which(x=="N")] <- NA
  return(y)
}

X <- apply(GBS[, -c(1:f.column)],1,parse.GBS)

setwd("/Users/ssapkot/Documents/Experiments/SAP_GS_PopStr/Results/Prediction_Results/WR_AR_SameTPsize/WR_Pred_Accu/")

for (a in 1:5){
    
    Y1 = Y[Y$Cluster == a,]
    nr = floor(nrow(Y1)/5)
    cvf <- data.frame(matrix("",nrow=nr*5, ncol=0))
    cluster <- Y1$Race[a]
    
    for (b in 1:100){
        
        c.list = vector("list", 5)
        
        ctemp = sample(Y1$SN,5*nr)
        c.list[[b]] = ctemp
        cvf[b] <- cbind(c.list[[b]])
    }
    fold <- rep(1:5, times = nr)
    print(head(cvf))
    
for (j in 5:12) {

  for (i in 1:101) {
  CV.fold <- paste("V",toString(i-1),sep='')
 
  if (CV.fold == "V0") {
    Total_Result <- c()
      result <- c()
  }
  else {  
    
  Z <- cvf[,CV.fold]
  Z <- sort(Z) ##sort randomly selected individuals by taxa order, and so the pheno and geno will be in the same order when subsetted
  
  X1 <- X[Z,]
  
  A <- A.mat(X1)
  
  rownames(A) <- 1:nrow(X1)
  P <- Y[Z,]
  cvs <- sample(fold)#fold is determined by folds assign while making the cvf dataframe 
  print(cvs)
  y = P[,j]
  col = names(P[j])
  
  yhat <- data.frame(cbind( y, yhat = 0))
  yhat$yhat <- as.numeric(yhat$yhat)
  row.names(yhat) <- row.names(y)
  
      corr <- c()
      var_x <- c()
      var_y <- c()
      cov_xy <- c()
  
      for (k in 1:5) {
    # Make training (TRN) and testing (TST) dfs
    tst <- which(cvs == k) ##cvs == whichever cluster/race is to be predicted
    yNA <- y
    yNA[tst] <- NA # Mask yields for validation set
    df <- data.frame(y=yNA,gid=1:nrow(A)) # Set up dataframe with traits and genotype labels (same order as in A1) 
    
    # Build rrBLUP model and save yhat for the masked values
    rrblup <- kin.blup(df,K=A,geno="gid",pheno="y") #optional parameters: fixed effects, gaussian kernel, covariates
    yhat$yhat[tst] = rrblup$pred[tst]
    }
    corr <- cor(yhat$y[tst],yhat$yhat[tst],use="complete")
    var_x <- var(yhat$yhat[tst], use="complete")
    var_y <- var(yhat$y[tst], use="complete")
    cov_xy <- cov(yhat$y[tst],yhat$yhat[tst], use="complete")
      
      result <- c(corr,var_x,var_y,cov_xy)
          
      }
      
      Total_Result <- cbind(Total_Result,result)
  
  }
    rownames(Total_Result) <- c("corr","var_x","var_y","cov_xy")
    
  write.csv(Total_Result, file = paste("WR_Corr-Cov_",a,"_",col,".csv", sep=""))
    }

}
