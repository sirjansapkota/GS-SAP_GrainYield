# Documentations for Sapkota _et. al._ 2019


## Scripts

#### CV1_prediction.R, CV2_AR_prediction.R, CV2_WR_prediction.R
Contains R Script to implement GBLUP and cross validations using _kin.blup_ function in R (rrBLUP).

#### Heterozygosity_Pegas.R
R script to calculate heterozygosity per site using R package _Pegas_.

####  LD_Calculation_HillandWeir.R
Calculations for expected values of r2 under drift equilibrium _Hill and Weir (1988)_. As implemented in _Remington, et al. (2001)_.

#### MCMCglmm_CV1_vCov.R
fit multi-response model in MCMCglmm to calculate variance-covariance components due to conditional expectations from race.

#### Stratified_Sampling_SAP_cluster.R
R script to create a cross validation file with individuals proportionally divided from each cluster into five equal folds.

#### vcf_to_matrix 
This script contains the functions _my.read.vcf_ to read SNP file in vcf format, and _parse.vcf_ function to create a genotype matrix (in -1,0,1 format) that can be used in the rrblup package.
