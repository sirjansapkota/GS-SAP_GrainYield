
rm(list=ls())
.libPaths(c('/Users/ssapkot/Documents/SAP/GenSel/MyRlibs',.libPaths()))

setwd("/Users/ssapkot/Documents/GenomicSelection1/Results/LD/")
ld <- read.table("SAP_mixed_LD_50_tassel_maf0.05.txt",stringsAsFactors = FALSE,header=TRUE) ## read in table with R^2 and distance calculated from TASSEL or wherever
ld_sub <- ld[ld$R.2 != "NaN",]
ld_sub$dist <- as.numeric(ld_sub$Dist_bp)
ld_sub2 <- ld_sub[ld_sub$dist != "NaN",]
ld_sub2$rsq <- ld_sub2$R.2

file <- ld_sub2[,c(1,2,7,8,15:19)]

#test if there are any NAs or NaN
#isTRUE(file$rsq=="NA")
#isTRUE(file$rsq=="NaN") 


## in my experience, changing the value of C does not do a whole lot of damage to your calculations
## you can try different values of Cs to see for yourself
## but reported C values range from about 0.5 to 2
Cstart <- c(C=0.1)

# let's fit a non linear model using the arbitrary C value, N is the number of the genotypes that have the SNP site
modelC <- nls(rsq ~ ((10+C*dist)/((2+C*dist)*(11+C*dist)))*(1+((3+C*dist)*(12+12*C*dist+(C*dist)^2))/(2*N*(2+C*dist)*(11+C*dist))), data=file, start=Cstart, control=nls.control(maxiter=100))

# extract rho, the recombination parameter in 4Nr
rho <- summary(modelC)$parameters[1]

# feed in the new value of rho to obtain LD values adjusted for their distances along the chromosome/genome
newrsq <- ((10+rho*file$dist)/((2+rho*file$dist)*(11+rho*file$dist)))*(1+((3+rho*file$dist)*(12+12*rho*file$dist+(rho*file$dist)^2))/(2*file$N*(2+rho*file$dist)*(11+rho*file$dist)))

newfile <- data.frame(file$dist, newrsq)

#maxld <- max(file$rsq) #using max LD value from initial input file
maxld <- max(newfile$newrsq,na.rm=TRUE) #using max LD value from adjusted data
halfdecay = maxld*0.5
halfdecaydist <- newfile$file.dist[which.min(abs(newfile$newrsq-halfdecay))]
newfile <- newfile[order(newfile$file.dist),]

# create the plot
plot(file$dist, file$rsq, pch=".", cex=2, xlab="distance (bp)", ylab="LD (r^2)")
lines(newfile$file.dist, newfile$newrsq, col="blue", lwd=2)
abline(h=0.2, col="red")
abline(v=halfdecaydist, col="green")
mtext(round(halfdecaydist,2), side=1, line=0.05, at=halfdecaydist, cex=0.75, col="green")

f1 <- data.frame(newfile$file.dist, newfile$newrsq)
xval <- f1[which.min(abs(0.2 - f1$newfile.newrsq)),] #find x value where y=0.2
xval[1,1]


library(ggplot2) #load the library

ggplot(rsq0.5dist10kb, aes(dist, rsq)) + geom_line()

###calculate mean and percentage for various statistics
mean(file$rsq, na.rm =T) ## mean of unadjusted rsquare
mean(newrsq,na.rm =T) ## mean of adjusted rsquared
mean(ld_sub$pDiseq,na.rm=T) ###mean of pvalues
mean(ld_sub$D.prime, na.rm = T) # mean of D'
r35 <- file[file$rsq > 0.35,] ##subset sites with r^2 > 0.35
100*nrow(r35)/nrow(file) ## % of sites with r^2 > 0.35
p4 <- ld_sub[ld_sub$pDiseq <0.0001,] ##subset sites with pvalue <0.0001
100*nrow(p4)/nrow(ld_sub) ## % of sites with pvalue <0.0001

########
