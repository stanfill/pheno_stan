library(rstan)
library(simplepheno)
library(data.table)
data(weather_temperate_2011);data(weather_heat_2011)
data(weather_temperate_2012);data(weather_heat_2012)

data(p_triangular)
data(phenology_temperate_2011);data(phenology_heat_2011)
data(phenology_temperate_2012);data(phenology_heat_2012)


tavg2011 <- (weather_temperate_2011$tmin+weather_temperate_2011$tmax)/2
havg2011 <- (weather_heat_2011$tmin+weather_heat_2011$tmax)/2
tmax2011 <- weather_temperate_2011$tmax
hmax2011 <- weather_heat_2011$tmax

tavg2012 <- (weather_temperate_2012$tmin+weather_temperate_2012$tmax)/2
havg2012 <- (weather_heat_2012$tmin+weather_heat_2012$tmax)/2
tmax2012 <- weather_temperate_2012$tmax
hmax2012 <- weather_heat_2012$tmax

temperate2011 <- phenology_temperate_2011
heat2011 <- phenology_heat_2011

temperate2012 <- phenology_temperate_2012
heat2012 <- phenology_heat_2012

tdoy2011 <- as.POSIXlt(weather_temperate_2011$date)$yday + 1
hdoy2011 <- as.POSIXlt(weather_heat_2011$date)$yday + 1

tdoy2012 <- as.POSIXlt(weather_temperate_2012$date)$yday + 1
hdoy2012 <- as.POSIXlt(weather_heat_2012$date)$yday + 1

####################
#Create a more informed data frame to analyze
temperate2011$year <- 2011
temperate2011$temp <- "Temperate"
heat2011$year <- 2011
heat2011$temp <- "Hot"
temperate2012$year <- 2012
temperate2012$temp <- "Temperate"
heat2012$year <- 2012
heat2012$temp <- "Hot"

allDat <- rbind(temperate2011,heat2011,temperate2012,heat2012)
allDat <- as.data.table(na.omit(allDat))
setkey(allDat,GID,year,temp)

redallDat <- allDat[GID!="GID6179253",]
redallDat <- redallDat[GID!="GID6179559",]
redallDat <- redallDat[GID!="GID5398160",]

redallDat <- redallDat[-seq(9,nrow(redallDat),by=9),]

table(redallDat$GID)
temp2011 <- redallDat[year==2011&temp=="Temperate"];setkey(temp2011,GID)
hot2011 <- redallDat[year==2011&temp=="Temperate"];setkey(hot2011,GID)
temp2012 <- redallDat[year==2012&temp=="Hot"];setkey(temp2012,GID)
hot2012 <- redallDat[year==2012&temp=="Hot"];setkey(hot2012,GID)


####################
#Triangular phenology model

#The first index is temp/year combination
#The second index is genotype (GID)
#The third index is observations for that temp-by-GID combination
ngids <- 27
gid_obs <- nrow(hot2011)/ngids

dthArray <- array(0,c(4,27,2))
dthArray[1,,] <- matrix(temp2011$dth,nrow=27)
dthArray[2,,] <- matrix(hot2011$dth,nrow=27)
dthArray[3,,] <- matrix(temp2012$dth,nrow=27)
dthArray[4,,] <- matrix(hot2012$dth,nrow=27)

dtmArray <- array(0,c(4,27,2))
dtmArray[1,,] <- matrix(temp2011$dtm,nrow=27)
dtmArray[2,,] <- matrix(hot2011$dtm,nrow=27)
dtmArray[3,,] <- matrix(temp2012$dtm,nrow=27)
dtmArray[4,,] <- matrix(hot2012$dtm,nrow=27)

#Each row corresponds to a year/temp combination
n_weatherObs <- 731
weatherMat <- matrix(c(tavg2011[1:n_weatherObs],havg2011[1:n_weatherObs],
                       tavg2012[1:n_weatherObs],havg2012[1:n_weatherObs]),nrow=4,byrow=TRUE)
doyMat <- matrix(c(tdoy2011[1:n_weatherObs],hdoy2011[1:n_weatherObs],
                   tdoy2012[1:n_weatherObs],hdoy2012[1:n_weatherObs]),nrow=4,byrow=TRUE)
tMaxMat <- matrix(c(tmax2011[1:n_weatherObs],hmax2011[1:n_weatherObs],
                    tmax2012[1:n_weatherObs],hmax2012[1:n_weatherObs]),nrow=4,byrow=TRUE)


pheno_dat_gid <- list(ndays=ncol(weatherMat), nobs=dim(dthArray)[3],ngid=dim(dthArray)[2],
                      nyears=nrow(weatherMat),
                      obs_tavg=weatherMat, doy = doyMat, obs_tmax=tMaxMat,
                      obs_dth=dthArray, obs_dtm=dtmArray,
                      tthLow=950, tthHigh=50,tthmLow=950, tthmHigh=50,
                      tlower=c(0,25,40), tupper=c(3,3,3))

initial_multigid <- function(){
  list(tmin=rnorm(1,0,1),topt=rnorm(1,25,1),tmax=rnorm(1,40,1),
       sigma_dth=runif(1,3,10),sigma_dtm=runif(1,3,10),
       tth_g=rnorm(dim(dthArray)[2],850,1),tthm_g=rnorm(dim(dthArray)[2],850,1),
       mu_tth=rnorm(1,850,1),sig_tth=runif(1,1,4),
       mu_tthm=rnorm(1,850,1),sig_tthm=runif(1,1,4))
}

multiGID_fit <- stan(file="multigid_pheno_tri.stan", data=pheno_dat_gid, algorithm="NUTS",
                     init=initial_multigid,iter=1000, chains=2)

multiGID_fit
plot(multiGID_fit)
traceplot(multiGID_fit)



multiGID_fit2 <- stan(fit=multiGID_fit,data=pheno_dat_gid,init=initial_multigid,iter=2500,chains=2)

multiGID_fit2
plot(multiGID_fit2)
traceplot(multiGID_fit2)
