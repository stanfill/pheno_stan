library(rstan)
load("~/PhilCollaboration/stan_results/TestResults_GID_Grouped.RData")

simDf <- as.data.frame(fit)
nDraws <- nrow(simDf)

estDth <- grep("dthHat",colnames(simDf))

residsMat <- matrix(0,nDraws,pheno_dat_gid$nobs)
obsMat <- data.frame(dth=pheno_dat_gid$obs_dth,year_t=pheno_dat_gid$year_temp,gidgp=pheno_dat_gid$gidgp)

for(i in 1:nDraws){
  dthEstMat <- matrix(data.matrix(simDf[i,estDth]),nrow=pheno_dat_gid$nyears)
  residsMat[i,] <- obsMat$dth- dthEstMat[cbind(obsMat$year_t,obsMat$gidgp)]
  
}

head(residsMat[,1:10])
