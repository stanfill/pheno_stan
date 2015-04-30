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
#Combine all of the year/temp. combinations
allDat <- rbind(temperate2011,heat2011,temperate2012,heat2012,temperate2013,heat2013)
allDat <- as.data.table(na.omit(allDat))
setkey(allDat,GID,year,temp)
#table(allDat$GID)

####################
#Use classification tree results to create GID groups
G1GID <- c("GID6056245", "GID6171893", "GID6174886", "GID6175172", "GID6176178", "GID6176346", "GID6176523", "GID6177599", "GID6179253", "GID6179559")
G2GID <- c("GID16122",   "GID2465",    "GID3895",    "GID5390612", "GID5398160", "GID5999777", "GID6000921", "GID6175024", "GID6178401", "GID6178783")
G3GID <- c("GID4556647", "GID5077000", "GID5325839", "GID5343246", "GID5397958", "GID5423688", "GID5995410", "GID6179128", "GID6179222", "GID775") 

allDat$GIDgp <- 1
allDat[GID%in%G2GID,]$GIDgp <- 2
allDat[GID%in%G3GID,]$GIDgp <- 3

ngids <- max(allDat$GIDgp)
####################
#Create year_temp groups
allDat$year_temp <- 1 #Hot 2011
allDat[year==2011&temp!="Hot"]$year_temp <- 2 #Temperate 2011
allDat[year==2012&temp=="Hot"]$year_temp <- 3 #Hot 2012
allDat[year==2012&temp!="Hot"]$year_temp <- 4 #Temperate 2012
allDat[year==2013&temp=="Hot"]$year_temp <- 5 #Hot 2013
allDat[year==2013&temp!="Hot"]$year_temp <- 6 #Temperate 2013

####################
#Create average temp, day of year (doy) and max temp matrix

#Each row corresponds to a year/temp combination
n_weatherObs <- 731
tavgMat <- matrix(c(havg2011[1:n_weatherObs],tavg2011[1:n_weatherObs],
                       havg2012[1:n_weatherObs],tavg2012[1:n_weatherObs],
                       havg2013[1:n_weatherObs],tavg2013[1:n_weatherObs]),nrow=6,byrow=TRUE)
doyMat <- matrix(c(hdoy2011[1:n_weatherObs],tdoy2011[1:n_weatherObs],
                   hdoy2012[1:n_weatherObs],tdoy2012[1:n_weatherObs],
                   hdoy2013[1:n_weatherObs],tdoy2013[1:n_weatherObs]),nrow=6,byrow=TRUE)
tMaxMat <- matrix(c(hmax2011[1:n_weatherObs],tmax2011[1:n_weatherObs],
                    hmax2012[1:n_weatherObs],tmax2012[1:n_weatherObs],
                    hmax2013[1:n_weatherObs],tmax2013[1:n_weatherObs]),nrow=6,byrow=TRUE)


pheno_dat_gid <- list(ndays=ncol(tavgMat), nobs=nrow(allDat),ngid=ngids,
                      nyears=nrow(tavgMat),
                      obs_tavg=tavgMat, doy = doyMat, obs_tmax=tMaxMat,
                      obs_dth=allDat$dth, year_temp=allDat$year_temp, gidgp=allDat$GIDgp,
                      tthLow=950, tthHigh=50,
                      tmin=0,tmax=45,topt=28)

initial_multigid <- function(){
  list(sigma_dth=runif(1,3,10),tth_g=rnorm(ngids,1000,25),ppsen=rep(70,ngids))
}


##########
#Run two chains in parallel

f1 <- stan(file="TempGIDSeparate/philset_dth_only_new_structure.stan",data=pheno_dat_gid,
           init=initial_multigid,chains=1, iter=1)
#seed <- 12345
num_core  <-  2
CL  <-  makeCluster(num_core, outfile = 'cluster.log')
clusterExport(cl = CL, c("pheno_dat_gid","initial_multigid", "f1","ngids"))
sflist1 <-parLapply(CL, 1:2,
               fun = function(i) {
                 require(rstan)
                 stan(fit = f1, data = pheno_dat_gid, init=initial_multigid,
                        chains = 1, chain_id = i, iter=500, refresh = -1)
                 })
fit <- sflist2stanfit(sflist1)
stopCluster(CL)
beep()

fit
traceplot(fit)

##########
#Run one chain at a time

multiGID_fit <- stan(file="TempGIDSeparate/philset_dth_only_new_structure.stan", data=pheno_dat_gid, algorithm="NUTS",
                     init=initial_multigid,iter=500, chains=2)

multiGID_fit <- stan(fit=multiGID_fit, data=pheno_dat_gid, algorithm="NUTS",
                     init=initial_multigid,iter=500, chains=2)

multiGID_fit
plot(multiGID_fit)
traceplot(multiGID_fit)
