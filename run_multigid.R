library(rstan)
library(simplepheno)
data(weather_temperate_2011);data(weather_heat_2011)
data(p_triangular)
data(phenology_temperate_2011);data(phenology_heat_2011)

tavg2011 <- (weather_temperate_2011$tmin+weather_temperate_2011$tmax)/2
havg2011 <- (weather_heat_2011$tmin+weather_heat_2011$tmax)/2
temperate2011 <- phenology_temperate_2011
heat2011 <- phenology_heat_2011

####################
#No genetic information include

#The first index is temperature
#The second index is genotype (GID)
#The third index is observations for that temp-by-GID combination
dthArray <- array(rep(0,120),c(2,2,30))
dthArray[1,,] <- matrix(temperate2011$dth,nrow=2)
dthArray[2,,] <- matrix(heat2011$dth,nrow=2)


weatherMat <- matrix(c(tavg2011[1:500],havg2011[1:500]),nrow=2,byrow=TRUE)

pheno_dat_gid <- list(ndays=ncol(weatherMat), nobs=dim(dthArray)[3],ngid=nrow(dthArray),nyears=2,
                      tavg=weatherMat, dths=dthArray,
                      tthLow=925, tthHigh=1375,
                      tlower=c(-5,20,30), tupper=c(5,30,50))

initial_multigid <- function(){
  list(theta1=rnorm(1,0,1),theta2=rnorm(1,25,1),theta3=rnorm(1,40,1),sigma=rnorm(1,3,1),
                    tthpar=rnorm(nrow(dthArray),950,1),mu_tth=rnorm(1,950,1),sig_tth=rnorm(1,1,1))
}

multiGID_fit <- stan(file="multigid_pheno_tri.stan",data=pheno_dat_gid,
                     init=initial_multigid,iter=5000,chains=2)

multiGID_fit
plot(multiGID_fit)
traceplot(multiGID_fit)
