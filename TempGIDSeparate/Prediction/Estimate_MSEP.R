###############
#Estimate squared bias based on Wallach et al. framework paper
library(rstan)
library(ggplot2)
library(reshape2)
load("~/PhilCollaboration/stan_results/Six_GID_Groups_One_PPSEN_Triangle_5000.RData")
triFit <- fit
load("~/PhilCollaboration/stan_results/Six_GID_Groups_One_PPSEN_Wang_5000.RData")
wangFit <- fit
load("~/PhilCollaboration/stan_results/Six_GID_Groups_One_PPSEN_Trapezoid_5000.RData")
trapFit <- fit


resultsMat <- data.frame(dth=pheno_dat_gid$obs_dth,year_t=pheno_dat_gid$year_temp,
                         gidgp=pheno_dat_gid$gidgp,GID=allDat$GID)
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
(meanTriError <- mean(biasTri))


#Take the results for the Wang model and estimate bias and model variance
simMatWang <- as.matrix(wangFit)
fHatXhatThetaHatWang <- simMatWang[,grep("dthHat",colnames(simMatWang))]
fHatXhatWang <- colMeans(fHatXhatThetaHatWang) #Average over parameter draws  
fHatXhatMatWang <- matrix(data.matrix(fHatXhatWang),nrow=pheno_dat_gid$nyears)
biasWang <- resultsMat$dth-fHatXhatMatWang[cbind(resultsMat$year_t,resultsMat$gidgp)]
resultsMat$Wang <- fHatXhatMatWang[cbind(resultsMat$year_t,resultsMat$gidgp)]
(sbHatWang <- mean(biasWang^2))
(modelVarWang <- var(as.numeric(fHatXhatThetaHatWang)))
(meanWangError <- mean(biasWang))


#Take the results for the Trapezoid model and estimate bias and model variance
simMatTrap <- as.matrix(trapFit)
fHatXhatThetaHatTrap <- simMatTrap[,grep("dthHat",colnames(simMatTrap))]
fHatXhatTrap <- colMeans(fHatXhatThetaHatTrap) #Average over parameter draws  
fHatXhatMatTrap <- matrix(data.matrix(fHatXhatTrap),nrow=pheno_dat_gid$nyears)
biasTrap <- resultsMat$dth-fHatXhatMatTrap[cbind(resultsMat$year_t,resultsMat$gidgp)]
resultsMat$Trap <- fHatXhatMatTrap[cbind(resultsMat$year_t,resultsMat$gidgp)]
(sbHatTrap <- mean(biasTrap^2))
(modelVarTrap <- var(as.numeric(fHatXhatThetaHatTrap)))
(meanTrapError <- mean(biasTrap))


##################
#Average over repeated measurements
library(plyr)
avgResultsMat <- ddply(resultsMat,.(GID,year_t),summarize,avgDTH=mean(dth),Triangle=mean(Triangle),
                       Wang=mean(Wang),Trap=mean(Trap))
avgTriError <- avgResultsMat$avgDTH-avgResultsMat$Triangle
mean(avgTriError)
mean(biasTri)
sd(avgTriError)
sd(biasTri)

avgWangError <- avgResultsMat$avgDTH-avgResultsMat$Wang
mean(avgWangError)
mean(biasWang)
sd(avgWangError)
sd(biasWang)

avgTrapError <- avgResultsMat$avgDTH-avgResultsMat$Trap
mean(avgTrapError)
mean(biasTrap)
sd(avgTrapError)
sd(biasTrap)

##################
#Organize the data to ease plotting

biasDF <- data.frame(Year=resultsMat$Year,Temp=resultsMat$Temp,GIDgp=resultsMat$gidgp)
biasDF$Triangle <- resultsMat$dth-resultsMat$Triangle
biasDF$Wang <- resultsMat$dth-resultsMat$Wang
biasDF$Trap <- resultsMat$dth-resultsMat$Trap

biasMelt <- melt(biasDF,id.vars=c("Year","Temp","GIDgp"))

qplot(variable,value,facets=Year~Temp,data=biasMelt,colour=factor(GIDgp),size=I(2))+
  geom_hline(yintercept=0,colour='gray50')+theme_bw()

################
#Decomponse model variance
#Values from Wallach framework paper
I <- 2 #number of models
J <- 6 #number of input vectors (year - sowing date combos)
K <- nDraws #Number of parameter vectors drawn
L <- length(pheno_dat_gid$obs_dth) #Target population (?)

library(lme4)
estsMat <- data.frame(Est=c(fHatXhatMatTri,fHatXhatMatWang,fHatXhatMatTrap),GID=rep(1:6,each=6),
                      Year=rep(c(2011,2012,2013),each=2),Temp=c("Hot","Temp"),
                      Model=rep(c("Tri","Wang","Trap"),each=18))
estsMat$GIDf <- as.factor(estsMat$GID)
estsMat$Yearf <- as.factor(estsMat$Year)
estsMat$Input <- paste(estsMat$Temp,estsMat$Yearf)


UQfit <- lm(Est~Input*Model,data=estsMat)
anova(UQfit)
summary(UQfit)

UQfit2 <- lm(Est~Yearf*Temp*Model,data=estsMat)
anova(UQfit2)
summary(UQfit2)

plot(fHatXhatMatTri,fHatXhatMatWang)

###########
#See if there is some pattern in error
library(ggplot2)
qplot(Triangle,dth,data=resultsMat,facets=year_t~gidgp)+geom_abline(yintercept=0,slope=1)
qplot(Wang,dth,data=resultsMat,facets=year_t~gidgp)+geom_abline(yintercept=0,slope=1)
qplot(Trap,dth,data=resultsMat,facets=year_t~gidgp)+geom_abline(yintercept=0,slope=1)



#####################################################################
#Old stuff
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
