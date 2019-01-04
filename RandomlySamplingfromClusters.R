setwd("/Users/ssapkot/Documents/GenomicSelection1/SAP_Yield/01_Data/")

  
  rm(list=ls())
  .libPaths(c('/Users/ssapkot/Documents/SAP/GenSel/MyRlibs',.libPaths()))
  
  Y <- read.csv('/Users/ssapkot/Documents/GenomicSelection1/SAP_Yield/01_Data/pheno_caudatum_all.csv', header=T)
  
  #### sample 38 individuals from each cluster and make a subset for across cluster prediction ##################
  
  cvf <- data.frame(matrix("",nrow=190, ncol=0))
  
    cv.list <- vector("list",5)
    
    for (j in 1:100){
    
    c.list = vector("list", 5)
    
    for (i in 1:5) {
      ctemp = Y[Y$Cluster==i,]
      total = nrow(ctemp)
      #if (!(is.integer(total))) { total = as.integer(total) +1 }
      x2 = sample(ctemp$TaxaID,38)
      #g = rep(seq(1,5), total)
      #g2 = g[1:length(x2)]
      #x = split(x2, g2)
      c.list[[i]] = x2
    }
      cvf[j] <- c(c.list[[1]],c.list[[2]],c.list[[3]],c.list[[4]],c.list[[5]])
    }
    cluster <- c(rep(1,38),rep(2,38),rep(3,38),rep(4,38),rep(5,38))
    cvf<- cbind(cluster,cvf)
    write.csv(cvf,"Sampling_across_population_proportionallySampled.csv")
###unlist the cv.list with randomly sampled genotypes  (n = 38 x 5) and put it in a data frame#####
    

  


#randomly sample 38 individuals from each subpopulation, and make a list of subpopulations sampled
#merge all four subpopulations sampled into a file
#assign each sampled subset in a column vector and make a matrix with different CV samples with rowname corresponding to CV fold number and column name corresponding to list of taxa in sampled group
