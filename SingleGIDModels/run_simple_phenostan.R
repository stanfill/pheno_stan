library(rstan)
library(simplepheno)
data(weather_temperate_2011)
data(p_triangular)
data(phenology_heat_2011)

wdat <- weather_temperate_2011
p2011 <- phenology_heat_2011

####################
#No genetic information included

pheno_dat <- list(ndays=length(wdat$tmin), nobs=nrow(p2011),
                  tavg=(wdat$tmin+wdat$tmax)/2, dths=p2011$dth,
                  settth=850,
                  tlower=c(-10,20,30), tupper=c(10,30,40))

initial <- function(){
  list(theta1=-5,theta2=30,theta3=35,sigma=3)
}

first_fit <- stan(file="Simple/first_pheno_tri.stan",data=pheno_dat,init=initial,iter=1000,chains=1)
first_fit
plot(first_fit)

####################
#Estimate tt.h parameter, one for all data, triangular phenology model

pheno_dat_gid <- list(ndays=length(wdat$tmin), nobs=nrow(p2011),
                  tavg=(wdat$tmin+wdat$tmax)/2, dths=p2011$dth,
                  tthLow=925, tthHigh=1375,
                  tlower=c(-5,20,30), tupper=c(5,30,50))

initial_gid <- function(){
  list(theta1=0,theta2=25,theta3=35,sigma=3,tthpar=950)
}

oneGID_fit <- stan(file="Simple/gid_pheno_tri.stan",data=pheno_dat_gid,init=initial_gid,iter=1000,chains=1)
oneGID_fit
plot(oneGID_fit)
traceplot(oneGID_fit)

####################
#Estimate tt.h parameter, one for all data, trapezoidal phenology model

pheno_dat_gid_trap <- list(ndays=length(wdat$tmin), nobs=nrow(p2011),
                        tavg=(wdat$tmin+wdat$tmax)/2, dths=p2011$dth,
                        tthLow=925, tthHigh=1375,
                        tlower=c(-5,20,20,30), tupper=c(5,30,30,50))

initial_gid_trap <- function(){
  list(theta1=0,theta2=23,theta3=27,theta4=35,sigma=3,tthpar=950)
}

oneGID_trap_fit <- stan(file="Simple/gid_pheno_trap.stan",data=pheno_dat_gid_trap,
                        init=initial_gid_trap,iter=1000,chains=1)
oneGID_trap_fit
plot(oneGID_trap_fit)
traceplot(oneGID_trap_fit)

####################
#Estimate tt.h parameter, one for all data, Wang Engel phenology model

pheno_dat_gid_wang <- list(ndays=length(wdat$tmin), nobs=nrow(p2011),
                           tavg=(wdat$tmin+wdat$tmax)/2, dths=p2011$dth,
                           tthLow=925, tthHigh=1375,
                           tlower=c(-5,20,30), tupper=c(5,30,50))

initial_gid_wang <- function(){
  list(theta1=0,theta2=23,theta3=35,sigma=3,tthpar=950)
}

oneGID_wang_fit <- stan(file="Simple/gid_pheno_wang.stan",data=pheno_dat_gid_wang,
                        init=initial_gid_wang,iter=1000,chains=1)

oneGID_wang_fit
plot(oneGID_wang_fit)
traceplot(oneGID_wang_fit)



