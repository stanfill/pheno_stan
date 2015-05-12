library(rstan)
load("~/PhilCollaboration/stan_results/Three_GID_Groups_One_PPSEN_Triangle_2500.RData")

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


################
I <- 2 #number of models
J <- 6 #number of input vectors (year - sowing date combos)
K <- nDraws #Number of parameter vectors drawn
L <- length(pheno_dat_gid$obs_dth) #Target population (?)

#Estimate squared bias based on Wallach et al. framework paper
library(rstan)
library(ggplot2)
library(reshape2)
load("~/PhilCollaboration/stan_results/Three_GID_Groups_One_PPSEN_Triangle_2500.RData")
triFit <- fit
load("~/PhilCollaboration/stan_results/Three_GID_Groups_One_PPSEN_Wang_2500.RData")
wangFit <- fit
load("~/PhilCollaboration/stan_results/Three_GID_Groups_One_PPSEN_Topt_Wang_2000.RData")
wangToptFit <- fit
load("~/PhilCollaboration/stan_results/Three_GID_Groups_Three_PPSEN_Triangle_2500.RData")
triThreePPFit <- fit

resultsMat <- data.frame(dth=pheno_dat_gid$obs_dth,year_t=pheno_dat_gid$year_temp,gidgp=pheno_dat_gid$gidgp)
resultsMat$Year <- 2011
resultsMat[resultsMat$year_t%in%c(3,4),]$Year <- 2012
resultsMat[resultsMat$year_t%in%c(5,6),]$Year <- 2013
resultsMat$Temp <- "Temperate"
resultsMat[resultsMat$year_t%in%c(1,3,5),]$Temp <- "Hot"

#Take the results for the Triangular model and estimate bias and model variance
simMatTri <- as.matrix(triFit)
fHatXhatThetaHatTri <- simMatTri[,grep("dthHat",colnames(simMatTri))]
fHatXhatTri <- colMeans(fHatXhatThetaHatTri) #Average over parameter draws  
fHatXhatMatTri <- matrix(data.matrix(fHatXhatTri),nrow=pheno_dat_gid$nyears)
biasTri <- resultsMat$dth-fHatXhatMatTri[cbind(resultsMat$year_t,resultsMat$gidgp)]
resultsMat$Triangle <- fHatXhatMatTri[cbind(resultsMat$year_t,resultsMat$gidgp)]
(sbHatTri <- mean(biasTri^2))
(modelVarTri <- var(as.numeric(fHatXhatThetaHatTri)))

#Take the results for the Triangular model with Three PPSEN and estimate bias and model variance
simMatTriThree <- as.matrix(triThreePPFit)
fHatXhatThetaHatTriThree <- simMatTriThree[,grep("dthHat",colnames(simMatTriThree))]
fHatXhatTriThree <- colMeans(fHatXhatThetaHatTriThree) #Average over parameter draws  
fHatXhatMatTriThree <- matrix(data.matrix(fHatXhatTriThree),nrow=pheno_dat_gid$nyears)
biasTriThree <- resultsMat$dth-fHatXhatMatTriThree[cbind(resultsMat$year_t,resultsMat$gidgp)]
resultsMat$TriThree <- fHatXhatMatTriThree[cbind(resultsMat$year_t,resultsMat$gidgp)]
(sbHatTriThree <- mean(biasTriThree^2))
(modelVarTriThree <- var(as.numeric(fHatXhatThetaHatTriThree)))

#Take the results for the Wang model (fixed topt) and estimate bias and model variance
simMatWang <- as.matrix(wangFit)
fHatXhatThetaHatWang <- simMatWang[,grep("dthHat",colnames(simMatWang))]
fHatXhatWang <- colMeans(fHatXhatThetaHatWang) #Average over parameter draws  
fHatXhatMatWang <- matrix(data.matrix(fHatXhatWang),nrow=pheno_dat_gid$nyears)
biasWang <- resultsMat$dth-fHatXhatMatWang[cbind(resultsMat$year_t,resultsMat$gidgp)]
resultsMat$Wang <- fHatXhatMatWang[cbind(resultsMat$year_t,resultsMat$gidgp)]
(sbHatWang <- mean(biasWang^2))
(modelVarWang <- var(as.numeric(fHatXhatThetaHatWang)))

#Take the results for the Wang model (estimated topt) and estimate bias and model variance
simMatWangTopt <- as.matrix(wangToptFit)
fHatXhatThetaHatWangTopt <- simMatWangTopt[,grep("dthHat",colnames(simMatWangTopt))]
fHatXhatWangTopt <- colMeans(fHatXhatThetaHatWangTopt) #Average over parameter draws  
fHatXhatMatWangTopt <- matrix(data.matrix(fHatXhatWangTopt),nrow=pheno_dat_gid$nyears)
biasWangTopt <- resultsMat$dth-fHatXhatMatWangTopt[cbind(resultsMat$year_t,resultsMat$gidgp)]
resultsMat$WangTopt <- fHatXhatMatWangTopt[cbind(resultsMat$year_t,resultsMat$gidgp)]
(sbHatWangTopt <- mean(biasWangTopt^2))
(modelVarWangTopt <- var(as.numeric(fHatXhatThetaHatWangTopt)))


##################
#Organize the data to ease plotting

biasDF <- data.frame(Year=resultsMat$Year,Temp=resultsMat$Temp,GIDgp=resultsMat$gidgp)
biasDF$Triangle <- resultsMat$dth-resultsMat$Triangle
biasDF$TriThree <- resultsMat$dth-resultsMat$TriThree
biasDF$Wang <- resultsMat$dth-resultsMat$Wang
biasDF$WangTopt <- resultsMat$dth-resultsMat$WangTopt

biasMelt <- melt(biasDF,id.vars=c("Year","Temp","GIDgp"))

qplot(variable,value,facets=Year~Temp,data=biasMelt,colour=factor(GIDgp),size=I(2))+
  geom_hline(yintercept=0,colour='gray50')+theme_bw()

################
#Decomponse model variance
library(lme4)
estsMat <- data.frame(Est=c(fHatXhatMatTri,fHatXhatMatWang),GID=rep(1:3,each=6),
                      Year=rep(c(2011,2012,2013),each=2),Temp=c("Hot","Temp"),
                      Model=rep(c("Tri","Wang"),each=18))
estsMat$GIDf <- as.factor(estsMat$GID)
estsMat$Yearf <- as.factor(estsMat$Year)
estsMat$Input <- paste(estsMat$Temp,estsMat$Yearf)

UQfit <- lm(Est~Input*Model,data=estsMat)
anova(UQfit)

plot(fHatXhatMatTri,fHatXhatMatWang)

###########
#See if there is some pattern in error
library(ggplot2)
qplot(Triangle,dth,data=resultsMat,facets=year_t~gidgp)+geom_abline(yintercept=0,slope=1)
qplot(Wang,dth,data=resultsMat,facets=year_t~gidgp)+geom_abline(yintercept=0,slope=1)
qplot(WangTopt,dth,data=resultsMat,facets=year_t~gidgp)+geom_abline(yintercept=0,slope=1)

