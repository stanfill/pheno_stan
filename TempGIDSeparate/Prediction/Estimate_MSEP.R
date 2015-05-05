library(rstan)
library(simplepheno)
load("Z:/bragg_pheno_stan/Results/Sep_GID_One_PPSEN_Topt_2500.RData")
source('TempGIDSeparate/Prediction/FastPhenoFun.R')

simDf <- as.data.frame(fit)
N_draws <- nrow(simMat)

tthcols <- grep("tth_g",colnames(simDf))

meani <- matrix(0,N_draws,ngids)

