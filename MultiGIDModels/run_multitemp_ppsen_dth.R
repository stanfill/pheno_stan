library(rstan)
library(simplepheno)
library(data.table)
library(parallel)
library(beepr)
data(weather_temperate_2011);data(weather_heat_2011)
data(weather_temperate_2012);data(weather_heat_2012)
data(weather_temperate_2013);data(weather_heat_2013)

data(p_triangular)
data(phenology_temperate_2011);data(phenology_heat_2011)
data(phenology_temperate_2012);data(phenology_heat_2012)
data(phenology_temperate_2013);data(phenology_heat_2013)

tavg2011 <- (weather_temperate_2011$tmin+weather_temperate_2011$tmax)/2
havg2011 <- (weather_heat_2011$tmin+weather_heat_2011$tmax)/2
tmax2011 <- weather_temperate_2011$tmax
hmax2011 <- weather_heat_2011$tmax

tavg2012 <- (weather_temperate_2012$tmin+weather_temperate_2012$tmax)/2
havg2012 <- (weather_heat_2012$tmin+weather_heat_2012$tmax)/2
tmax2012 <- weather_temperate_2012$tmax
hmax2012 <- weather_heat_2012$tmax

tavg2013 <- (weather_temperate_2013$tmin+weather_temperate_2013$tmax)/2
havg2013 <- (weather_heat_2013$tmin+weather_heat_2013$tmax)/2
tmax2013 <- weather_temperate_2013$tmax
hmax2013 <- weather_heat_2013$tmax


temperate2011 <- phenology_temperate_2011
temperate2011$year <- 2011
temperate2011$temp <- "Temperate"
heat2011 <- phenology_heat_2011
heat2011$year <- 2011
heat2011$temp <- "Hot"
tdoy2011 <- as.POSIXlt(weather_temperate_2011$date)$yday + 1
hdoy2011 <- as.POSIXlt(weather_heat_2011$date)$yday + 1


temperate2012 <- phenology_temperate_2012
temperate2012$year <- 2012
temperate2012$temp <- "Temperate"
heat2012 <- phenology_heat_2012
heat2012$year <- 2012
heat2012$temp <- "Hot"
tdoy2012 <- as.POSIXlt(weather_temperate_2012$date)$yday + 1
hdoy2012 <- as.POSIXlt(weather_heat_2012$date)$yday + 1


temperate2013 <- phenology_temperate_2013
temperate2013$year <- 2013
temperate2013$temp <- "Temperate"
heat2013 <- phenology_heat_2013
heat2013$year <- 2013
heat2013$temp <- "Hot"
tdoy2013 <- as.POSIXlt(weather_temperate_2013$date)$yday + 1
hdoy2013 <- as.POSIXlt(weather_heat_2013$date)$yday + 1


####################
#Create a more informed data frame to analyze
allDat <- rbind(temperate2011,heat2011,temperate2012,heat2012,temperate2013,heat2013)
allDat <- as.data.table(na.omit(allDat))
setkey(allDat,GID,year,temp)
#table(allDat$GID)

redallDat <- allDat[GID!="GID6179253",]
redallDat <- redallDat[GID!="GID6179559",]
redallDat <- redallDat[GID!="GID5398160",]

redallDat <- redallDat[-seq(14,nrow(redallDat),by=14),]
redallDat <- redallDat[-seq(9,nrow(redallDat),by=13),]
#table(redallDat$GID)

temp2011 <- redallDat[year==2011&temp=="Temperate"];setkey(temp2011,GID)
hot2011 <- redallDat[year==2011&temp=="Hot"];setkey(hot2011,GID)
temp2012 <- redallDat[year==2012&temp=="Temperate"];setkey(temp2012,GID)
hot2012 <- redallDat[year==2012&temp=="Hot"];setkey(hot2012,GID)
temp2013 <- redallDat[year==2013&temp=="Temperate"];setkey(temp2013,GID)
hot2013 <- redallDat[year==2013&temp=="Hot"];setkey(hot2013,GID)

####################
#Triangular phenology model

#The first index is temp/year combination
#The second index is genotype (GID)
#The third index is observations for that temp-by-GID combination
gid_obs <- nrow(hot2011)
nyears <- 6
dthArray <- matrix(0,6,gid_obs)
dthArray[1,] <- temp2011$dth
dthArray[2,] <- hot2011$dth
dthArray[3,] <- temp2012$dth
dthArray[4,] <- hot2012$dth
dthArray[5,] <- temp2013$dth
dthArray[6,] <- hot2013$dth


#Each row corresponds to a year/temp combination
n_weatherObs <- 731
weatherMat <- matrix(c(tavg2011[1:n_weatherObs],havg2011[1:n_weatherObs],
                       tavg2012[1:n_weatherObs],havg2012[1:n_weatherObs],
                       tavg2013[1:n_weatherObs],havg2013[1:n_weatherObs]),nrow=6,byrow=TRUE)
doyMat <- matrix(c(tdoy2011[1:n_weatherObs],hdoy2011[1:n_weatherObs],
                   tdoy2012[1:n_weatherObs],hdoy2012[1:n_weatherObs],
                   tdoy2013[1:n_weatherObs],hdoy2013[1:n_weatherObs]),nrow=6,byrow=TRUE)
tMaxMat <- matrix(c(tmax2011[1:n_weatherObs],hmax2011[1:n_weatherObs],
                    tmax2012[1:n_weatherObs],hmax2012[1:n_weatherObs],
                    tmax2013[1:n_weatherObs],hmax2013[1:n_weatherObs]),nrow=6,byrow=TRUE)


pheno_dat_gid <- list(ndays=ncol(weatherMat), nobs=gid_obs,
                      nyears=nrow(weatherMat),
                      obs_tavg=weatherMat, doy = doyMat, obs_tmax=tMaxMat,
                      obs_dth=dthArray, 
                      tthLow=950, tthHigh=50,
                      tlower=c(0,25,40), tupper=c(3,3,3))

initial_multigid <- function(){
  list(tmin=runif(1,-5,5),topt=runif(1,20,30),tmax=runif(1,30,50),
       sigma_dth=runif(1,3,10),
       tth_g=rnorm(nyears,1000,10),
       #mu_tth=rnorm(1,1000,10),sig_tth=runif(1,1,4),
       ppsen=rnorm(1,30,5))
}


##########
#Run two chains in parallel

f1 <- stan(file="multitemp_wt_ppsen_dth.stan",data=pheno_dat_gid,init=initial_multigid,chains=1, iter=1)
#seed <- 12345
num_core  <-  2
CL  <-  makeCluster(num_core, outfile = 'cluster.log')
clusterExport(cl = CL, c("pheno_dat_gid","initial_multigid", "f1","nyears"))
sflist1 <-parLapply(CL, 1:2,
               fun = function(i) {
                 require(rstan)
                 stan(fit = f1, data = pheno_dat_gid, init=initial_multigid,
                        chains = 1, chain_id = i, iter=1000, refresh = -1)
                 })
fit <- sflist2stanfit(sflist1)
stopCluster(CL)
beep()

fit
traceplot(fit)

##########
#Run one chain at a time

multiGID_fit <- stan(file="multitemp_wt_ppsen_dth.stan", data=pheno_dat_gid, algorithm="NUTS",
                     init=initial_multigid,iter=1000, chains=2)

multiGID_fit
plot(multiGID_fit)
traceplot(multiGID_fit)
