# Script to decompose CV1 covariance into expectation over races and covariance from individuals within race.

library(MCMCglmm)

# 'x' is  predicted value  and 'y' is observed phenotypic value for a given 
# trait, and 'race' represents the racial cluster to which the accessions belong

mTrait <- MCMCglmm(cbind(y,x) ~ trait - 1, random = ~idh(trait):race, 
                   rcov=~us(trait):units, data = data, 
                   family=c("gaussian","gaussian"), verbose = FALSE)

# total covariance calculation
cov_xy = cov(data$y.pred, data$y.obs)

# variance-covariance components due to individuals within a race
test  <- summary(m_DTA)$Rcovariances

# Cov(x, y) = Erace [Cov(x, y|race)] + Covrace [E(x|race), E(y|race)]
# expectation due to race
E_race <- cov_xy - test[2]
