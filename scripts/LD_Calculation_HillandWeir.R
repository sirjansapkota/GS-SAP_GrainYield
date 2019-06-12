# Calculations for expected values of r2 under drift equilibrium [Hill and Weir (1988)
# implemented in Remington, et al. (2001)]

# Calculate and save LD statistics from TASSEL using sliding window of 50 SNPs and maf = 0.1
# read in table with R^2 and distance calculated from TASSEL or wherever
ld <- read.table("SAP_mixed_LD_50_tassel_maf0.1.txt",stringsAsFactors = FALSE,header=TRUE) 
##remove sites that have missing values for distance or r2
ld_sub <- ld[ld$R.2 != "NaN",]
ld_sub$dist <- as.numeric(ld_sub$Dist_bp)
ld_sub2 <- ld_sub[ld_sub$dist != "NaN",]
ld_sub2$rsq <- ld_sub2$R.2

file <- ld_sub2[,c(1,2,7,8,15:19)]


# C values range from about 0.5 to 2, start with 0.1
Cstart <- c(C=0.1)

# fit a non linear model using the arbitrary C value, 
# N is the number of the genotypes that have the SNP site
modelC <- nls(rsq ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
                ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( 2*N*(2+C*dist) * (11+C*dist) ) ), 
              data=file, start=Cstart, control=nls.control(maxiter=100))

# extract rho, the recombination parameter in 4Nr
rho <- summary(modelC)$parameters[1]

# feed in the new value of rho to obtain LD values adjusted for their distances along the chromosome/genome
newrsq <- ( (10+rho*file$dist) / ( (2+rho*file$dist) * (11+rho*file$dist) ) ) *
  ( 1 + ( (3+rho * file$dist) * (12+12*rho*file$dist + (rho*file$dist)^2) ) / 
      (2*file$N*(2+rho*file$dist) * (11+rho*file$dist) ) )

newfile <- data.frame(file$dist, newrsq)

maxld <- max(newfile$newrsq,na.rm=TRUE) #using max LD value from adjusted data
halfdecay = maxld*0.5
halfdecaydist <- newfile$file.dist[which.min(abs(newfile$newrsq-halfdecay))]
newfile <- newfile[order(newfile$file.dist),]

# plotting the values
plot(file$dist, file$rsq, pch=".", cex=2, xlab="distance (bp)", ylab="LD (r^2)")
lines(newfile$file.dist, newfile$newrsq, col="blue", lwd=2)
abline(h=0.2, col="red")
abline(v=halfdecaydist, col="green")
mtext(round(halfdecaydist,2), side=1, line=0.05, at=halfdecaydist, cex=0.75, col="green")
