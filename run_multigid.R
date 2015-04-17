library(rstan)
library(simplepheno)
data(weather_temperate_2011);data(weather_heat_2011)
data(weather_temperate_2012);data(weather_heat_2012)

data(p_triangular)
data(phenology_temperate_2011);data(phenology_heat_2011)

tavg2011 <- (weather_temperate_2011$tmin+weather_temperate_2011$tmax)/2
havg2011 <- (weather_heat_2011$tmin+weather_heat_2011$tmax)/2
temperate2011 <- phenology_temperate_2011
heat2011 <- phenology_heat_2011

tdoy2011 <- as.POSIXlt(weather_temperate_2011$date)$yday + 1
hdoy2011 <- as.POSIXlt(weather_heat_2011$date)$yday + 1


####################
#Triangular phenology model

#The first index is temperature
#The second index is genotype (GID)
#The third index is observations for that temp-by-GID combination
dthArray <- array(rep(0,120),c(2,2,30))
dthArray[1,,] <- matrix(temperate2011$dth,nrow=2)
dthArray[2,,] <- matrix(heat2011$dth,nrow=2)

dtmArray <- array(rep(0,120),c(2,2,30))
dtmArray[1,,] <- matrix(temperate2011$dtm,nrow=2)
dtmArray[2,,] <- matrix(heat2011$dtm,nrow=2)

weatherMat <- matrix(c(tavg2011[1:500],havg2011[1:500]),nrow=2,byrow=TRUE)

doyMat <- matrix(c(tdoy2011[1:500],hdoy2011[1:500]),nrow=2,byrow=TRUE)

pheno_dat_gid <- list(ndays=ncol(weatherMat), nobs=dim(dthArray)[3],ngid=nrow(dthArray),nyears=2,
                      tavg=weatherMat, doy = doyMat, obs_dth=dthArray, obs_dtm=dtmArray,
                      tthLow=925, tthHigh=1375,tthmLow=925, tthmHigh=1375,
                      tlower=c(-5,20,30), tupper=c(5,30,50))

initial_multigid <- function(){
  list(tmin=rnorm(1,0,1),topt=rnorm(1,25,1),tmax=rnorm(1,40,1),
       sigma_dth=rnorm(1,3,1),sigma_dtm=rnorm(1,3,1),
       tth_g=rnorm(nrow(dthArray),950,1),tthm_g=rnorm(nrow(dthArray),950,1),
       mu_tth=rnorm(1,950,1),sig_tth=abs(rnorm(1,3,1)),mu_tthm=rnorm(1,950,1),sig_tthm=abs(rnorm(1,3,1)))
}

multiGID_fit <- stan(file="multigid_pheno_wang.stan",data=pheno_dat_gid,
                     init=initial_multigid,iter=1000,chains=1)

multiGID_fit
plot(multiGID_fit)
traceplot(multiGID_fit)


####################
#Wang-engel phenology model

multiGID_fit_wang <- stan(file="multigid_pheno_wang.stan",data=pheno_dat_gid,
                     init=initial_multigid,iter=1000,chains=1)

multiGID_fit_wang2 <- stan(fit=multiGID_fit_wang,data=pheno_dat_gid,
     init=initial_multigid,iter=2500,chains=5)

multiGID_fit_wang2
plot(multiGID_fit_wang2)
traceplot(multiGID_fit_wang2)
