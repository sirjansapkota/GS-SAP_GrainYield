
.libPaths(c('/Users/ssapkot/Documents/Experiments/MyRlibs',.libPaths()))
library(BGLR)

load("/Users/ssapkot/Documents/Experiments/SAP_GS_PopStr/VarianceComp/VarComp_SAP_GS_heritability_WPvsBP.RData")
load('/Users/ssapkot/Documents/Experiments/SAP_GS_PopStr/VarianceComp/programs/RKHS.rda') ## this contain RKHS function used to calculate variance for each eigenvectors and also eigenvalues

for (i in 1:ncol(Y)) {
setwd("/Users/ssapkot/Documents/Experiments/SAP_GS_PopStr/VarianceComp/output/") 

trait <- colnames(Y)    
dir.create(paste(trait[i]))
setwd(paste(trait[i]))

nIter <- 37000 #no. of iterations
burnIn <- 2000 # no. of burIn, outcomes of runs from burnIn are discarded

## Eigen-value decomposition of Genomic Relationship Matrix
EVD<-eigen(G)

for (nPCs in 0:20) { # run for loop to calculate within and between subpopulation variance each each eigenvector upto 20
  
dir.create(paste(nPCs,'_PCs',sep=''))
setwd(paste(nPCs,'_PCs','/',sep=''))

y<-Y[,i] ## Do this for each phenotypic values to calculated decomposition by traits

if(nPCs>0){
  PCs<-as.matrix(EVD$vectors[,1:nPCs])
  fm<-RKHS(y=y,XF=PCs,K=list(list(V=EVD$vectors,d=EVD$values,df0=5,S0=3/2)),nIter=nIter,burnIn=burnIn)
}else{
  fm<-RKHS(y=y,K=list(list(V=EVD$vectors,d=EVD$values,df0=5,S0=3/2)),nIter=nIter,burnIn=burnIn)
}
save(fm,file='fm.rda')
rm(fm)
    
setwd('../')
}
    
}

load("/Users/ssapkot/Documents/Experiments/SAP_GS_PopStr/VarianceComp/VarComp_SAP_GS_heritability_WPvsBP.RData")
for (i in 1:ncol(Y)) {
setwd("/Users/ssapkot/Documents/Experiments/SAP_GS_PopStr/VarianceComp/output/Using_G") 

trait <- colnames(Y)
setwd(paste(trait[i]))
    
nPCs<-0:20
VARE<-rep(NA,length(nPCs))
names(VARE)<-paste('PC',nPCs,SEP='')

VARU<-VARE
H2<-VARE
R2<-VARE
VARU_B<-cbind(VARE[-1],VARE[-1])
colnames(VARU_B)<-c('withPCs_B','withoutPCS_B')
VARU_W<-VARU_B
colnames(VARU_W)<-c('withPCs_W','withoutPCs_W')

row<-1
    
for(j in nPCs){   
  setwd(paste(j,'_PCs','/',sep=''))
  load('fm.rda')
  VARE[(1+j)]<-fm$varE
  VARU[(1+j)]<-fm$K[[1]]$varU
  R2[(j+1)]<-cor(Y[,1],fm$yHat)^2
  if(j==0){
     for(k in 1:20){
         VARU_B[k,2]<-fm$K[[1]]$cumMSa[k]
         VARU_W[k,2]<-(max(fm$K[[1]]$cumMSa)-fm$K[[1]]$cumMSa[k])
     } 
  }else{
    VARU_B[row,1]<-fm$K[[1]]$cumMSa[(j+1)]
    VARU_W[row,1]<-(max(fm$K[[1]]$cumMSa)-fm$K[[1]]$cumMSa[(j+1)])
    row<-row+1
  }
  setwd('../')

}

H2<-VARU/(VARE+VARU)
H2W<-VARU_W
H2W[]<-NA
H2W[,2]<-VARU_W[,2]/(VARU_W[,2]+VARE[1])
H2W[,1]<-VARU_W[,1]/(VARU_W[,1]+VARE[-1])

h20<-VARU[1]/(VARU[1]+VARE[1])
H2W<-rbind(rep(h20,2),H2W)


pdf('figures_EigenDecomp.pdf')
# Scatterplot showing PC1 against PC2
plot(EVD$vectors[,2:3],xlab='Eiegenvector 2',ylab='Eigenvector 1',main='PC1 vs PC2')

# Figure showing cumulative variation explained by number of eigenvalues
plot(c(0,cumsum(EVD$values))/389,x=0:389,xlab='Number of eigenvalues',ylim=c(0,1),type='l')
abline(a=0,b=1/389,lty=2) ; abline(h=1,lty=2)
dev.off()

    
pdf('figures.pdf')
 # Figure showing posterior means of within-subpopulation GenHeritability with effects of PCs
   tmp<-range(as.vector(H2W))
   plot(H2W[,2]~I(0:20),pch=15,col='red',xlab='Number of eigenvectors',ylab='Within-group heritability',type='o',ylim=tmp*c(.95,1.05))
   abline(h=range(H2W[,2]),lty=2)
   lines(x=0:20,y=H2W[,1],lty=1,col=4)
   points(x=0:20,y=H2W[,1],pch=15,col=4)
   abline(v=0)
   # labels
   lines(x=c(2,3),y=c(1,1),col='red')
   points(y=1,x=3.5,pch=15,col='red')
   text(x=5,cex=.9,y=1,col='red', label='Without Eigen Vectors as Fixed Effects')

   lines(x=c(2,3),y=c(.98,.98),col='blue')
   points(y=.98,x=3.5,pch=15,col='blue')
   text(x=5,cex=.9,y=.98,col='blue',label='With Eigen Vectors as Fixed Effects')
 # MCMC estimates of Posterior probabilites for each trait
   load('0_PCs/fm.rda')
   plot(fm$K[[1]]$probH1,type='l',col=4,ylab='Posterior Prob',xlab='Eigenvector')
   points(x=1:389,y=fm$K[[1]]$probH1,col=2,cex=.6)
dev.off()


VarE = as.data.frame(VARE[2:21])

Variances = as.data.frame(cbind(VARU_B[,1],VARU_W[,1], VarE[,1]))
colnames(Variances) = c("VarB","VarW","VarE")

    h2gA = c()
    h2gW = c()
    h2g  = c()
for (i in 1:nrow(Variances)){
  total  = sum(Variances[i,1]+Variances[i,2]+Variances[i,3])
  
  h2gA[i] = Variances[i,1]/total
  h2gW[i] = Variances[i,2]/total
  h2g[i] = h2gA[i] + h2gW[i]

}
VarHerit  =  cbind(Variances,h2gA,h2gW,h2g)
write.csv(VarHerit, file = paste("HeritabilityAndVarianceComp",trait[i],".csv",sep = "_"))
    
print(paste(trait[i],"_h2gA = ", mean(h2gA),":",sd(h2gA)),sep="")
print(paste(trait[i],"_h2gW = ", mean(h2gW),":",sd(h2gW)),sep="")
print(paste(trait[i],"_h2g = ", mean(h2g),":",sd(h2g)),sep="")
print(paste(trait[i],"_proportion_h2gA = ", mean(h2gA)/mean(h2g)),sep = "")
print(paste(trait[i],"_proportion_h2gW = ", mean(h2gW)/mean(h2g)),sep="")
    
    
}
