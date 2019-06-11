# Script to implement GBLUP and cross validations using kin.blup (rrBLUP)
## Cross validation was done using stratified/proportional sampling of accessions from different subpopulations/races in the diversity panel(SAP)
## read in the genotypic matrix 'X' and trait BLUEs 'Y' and 'cvf' file for cross validation folds

library(rrBLUP)

A = A.mat(X) #calculate additive/genomic relationship matrix
rownames(A) <- 1:nrow(X)

Pred_accuracy <- c()
Total_predValues <- c()

for (j in 1:ncol(Y)) {
corr <- c()
    PV = c()
  for (i in 1:101) {
    CV.fold <- paste("CV_",toString(i-1),sep='')
    
  if (CV.fold == "CV_0") {
    Pred_values <- c()
    corr <- c()
    }
  else {
    cvs <- cvf[,CV.fold]
  y <- Y[,j] #get vector of phenotypic values for jth phenotype
    colname = colnames(Y[j])
      
  yhat <- data.frame(cbind( y, yhat = 0))
  yhat$yhat <- as.numeric(yhat$yhat)
  row.names(yhat) <- row.names(y)
  accuracy <- c() 
  
  for(k in 1:5){ 
    # Make training (TRN) and testing (TST) dfs
    tst <- which(cvs == k)
    yNA <- y
    yNA[tst] <- NA # Mask yields for validation set
    df <- data.frame(y=yNA,gid=1:nrow(A)) # Set up dataframe with traits and genotype labels
    # Build rrBLUP model and save yhat for the masked values
    rrblup <- kin.blup(df,K=A,geno="gid",pheno="y")
     yhat$yhat[tst] <- rrblup$pred[tst]
  }
  accuracy <- cor(yhat$y, yhat$yhat, use = "complete.obs")
    Pred_values <- cbind(Pred_values,yhat$yhat)
      corr <-  rbind(corr, accuracy)
  }
      PV = cbind(Y$Cluster, Pred_values)
}

Pred_accuracy[[j]] <- corr
write.csv(PV,file=paste("Predicted_Values",colname,"csv",sep="."), row.names=rowname)
Total_predValues[[j]] <- Pred_values

}
write.csv(Pred_accuracy,"Prediction_accuracy_CV1.csv")


