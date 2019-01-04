library(rrBLUP)
A <- A.mat(X) ## Compute additive relationship matrix
rownames(A) <- 1:nrow(X)

Y <- read.csv('/Users/ssapkot/Documents/GenomicSelection1/SAP_Yield/01_Data/pheno_all.csv', header=T) #phenotype dataframe

cvs <- read.csv('/Users/ssapkot/Documents/GenomicSelection1/SAP_Yield/01_Data/CVFs_all.csv', row.names=1) ##cross validation file for sampling

  accuracy <- c()
  CV_accuracy <- c()
  
  y = Y$PH # read single atomic vector of phenotypic data from phenotype data frame
  
  for (j in 1:31) {
    CV.fold= paste('cv_', toString(j-1), sep='')
    df_out <- c()
    
    if(CV.fold =='cv_0'){
      #fit model to the entire data set and save model
      df <- data.frame(y=y,gid=1:nrow(A)) # Set up dataframe with traits and genotype labels (same order as in A1) 
      rrblup <- kin.blup(df,K=A,geno= "gid",pheno="y") #optional parameters: fixed effects, gaussian kernel, covariates
      save(rrblup,file='PH_all_fullmodel.RData')
    }
    
    else{
      # fit the model using the training sets designated by CV.
      tst=cvs[,CV.fold]
      yhat <- data.frame(cbind(y,yt, yhat = 0)) #make dataframe to with observed values and to store predicted values
      yhat$yhat <- as.numeric(yhat$yhat)
      row.names(yhat) <- row.names(y)
      
      for(k in 1:5){
        # Make training (TRN) and testing (TST) dfs
        test <- which(tst==k)
        yNA <- y
        yNA[test] <- NA # Mask yields for validation set
        df <- data.frame(y=yNA,gid=1:nrow(A)) # Set up dataframe with traits and genotype labels (same order as in A1) 
        
        # Build rrBLUP model and save yhat for the masked values
        rrblup <- kin.blup(df,K=A,geno="gid",pheno="y") #optional parameters: fixed effects, gaussian kernel, covariates
        yhat$yhat[test] <- rrblup$g[test]
        
      }
      accuracy <- cor(yhat$y, yhat$yhat) #prediction accuracy for one run of cross validation
      
      df_out <- data.frame(accuracy)
      
    }
    CV_accuracy <- rbind(CV_accuracy,df_out) #store all 30 iterations of prediction accuracies in a vectorized data frame
  }
