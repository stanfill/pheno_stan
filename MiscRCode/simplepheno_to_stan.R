#stan needs everything in a loop so I will turn his functions into a loop

tbase <- t.base <- 0
topt <- t.opt <- 25
tmax <- t.max <- 40
ppsen <- 100
ppthr <- 20
temp <- rnorm(10,27,1)


step <- function(x){
  ifelse(x>0,return(1),return(0))
}

pow <- function(x,y){
  return(x^y)
}

if_else <- function(test,yes,no){
  return(ifelse(test,yes,no))
}

cumulative_sum <- function(x){
  return(cumsum(x))
}

############### 
#His triangular function

triangular <- function(temp,base,opt,max)
  {
    val = vector(length=length(temp))
    
    val[temp<=opt] = (temp[temp<=opt] - base) / (opt - base)
    val[temp>opt] = (max - temp[temp>opt]) / (max-opt)
    val[val>1] = 1
    val[val<0] = 0
    
    return(val)
  }

#stan Triangular function
stan_triangular <- function(temp,t.base,t.opt,t.max){
  N <- length(temp)
  val <- z <- rep(0,N)
  for(i in 1:N){
    z[i] <- 1-step(temp[i]-t.opt) #Zi=1 if temp[i]<=t.opt; 0 otherwise
    val[i] <- max(min( z[i]*(temp[i]-t.base)/(t.opt-t.base)+(1-z[i])*(t.max-temp[i])/(t.max-t.opt),1),0) 
  }
  return(val)
}

triangular(temp,t.base,t.opt,t.max)-stan_triangular(temp,t.base,t.opt,t.max) #They match

############### 
#His trapezoidal function
temp <- rnorm(10,27,1)
opt1 <- 23
opt2 <- 27
trapezoidal <-function (temp, base, opt1, opt2, max) 
  {
    val = vector(length=length(temp))
    
    val[temp<=opt1] = (temp[temp<=opt1] - base) / (opt1 - base)
    val[temp>opt2] = (max - temp[temp>opt2]) / (max-opt2)
    val[val>1] = 1
    val[val<0] = 0
    
    return(val)
  }


#My trapezoidal function
N <- length(temp)
valp <- zp <- rep(NA,N)
for(i in 1:N){
  
  valp[i] <- if_else(temp[i]<=opt1,(temp[i] - t.base) / (opt1 - t.base),0)
  valp[i] <- if_else(temp[i]>opt2,(t.max - temp[i]) / (t.max-opt2),0)
  valp[i] <- max(min(valp[i],1),0)
}

trapezoidal(temp,t.base,opt1,opt2,t.max)-valp #They match

##################
wang_engel <-function(temp,base,opt,max){
    val = vector(length=length(temp))
    
    a = log(2) / log((max - base) / (opt - base))
    val = (2*(temp - base)^{a}*(opt - base)^{a} - (temp - base)^{2*a}) /(opt - base)^{2*a}
    val[temp<base|temp>max] = 0
    
    return(val)
}

stan_wang <-function(temp,tbase,topt,tmax){
    N <- length(temp)
    val = rep(0,N)
    a = log(2) / log((tmax - tbase) / (topt - tbase))
    denom <- pow(topt - tbase,2*a)
    
    for(i in 1:N){
      #val[i] = (2*(temp[i] - base)^{a}*(opt - base)^{a} - (temp[i] - base)^{2*a}) /(opt - base)^{2*a}
      
      val[i] = (2*pow(temp[i] - tbase,a)*pow(topt - tbase,a) - pow(temp[i] - tbase,2*a)) / denom
      
      val[i] <- if_else(temp[i]<tbase||temp[i]>tmax,0,val[i])
    }
    return(val)
}

wang_engel(temp,t.base,t.opt,t.max)-stan_wang(temp,tbase,topt,tmax)

#####################
#Daylength adjustment factor

id <- weather_temperate_2011
p <- list(ppsen=70)

##His function
doy = as.POSIXlt(id$date)$yday + 1

ppsen_fun <- function(doy,lat,ppsen){
  
  s1 = sin(lat * 0.01745)
  c1 = cos(lat * 0.01745)
  dec = 0.4093 * sin(0.0172 * (doy - 82.2))
  dlv = ((-s1 * sin(dec) - 0.1047)/(c1 * cos(dec)))
  dlv[dlv < -0.87] = -0.87
  dayl = 7.639 * acos(dlv)
  dayl.fac = 1 - ppsen/10000 * (20 - dayl)^2
  
  return(dayl.fac)
}
#My stan friendly version
ppsen <- 70
doy  <-  as.POSIXlt(id$date)$yday + 1
dayl_fac <- rep(0,length(doy))
lat <- id$lat
  
for(i in 1:length(doy)){
  s1 = sin(lat * 0.01745)
  c1 = cos(lat * 0.01745)
  dec = 0.4093 * sin(0.0172 * (doy[i] - 82.2))
  dlv = ((-s1 * sin(dec) - 0.1047)/(c1 * cos(dec)))
  dlv <- if_else(dlv< -0.87,-0.87,dlv)

  dayl = 7.639 * acos(dlv)
  dayl_fac[i] = 1 - ppsen/10000 * (20 - dayl)^2
  
}

#####################
#Vernelasation sensitivity adjustment

id <- weather_temperate_2011
tavg <- (id$tmin+id$tmax)/2
pbase <- -5; popt <- 7; pmax <- 15
tmax <- id$tmax

final_vrn_fac  <-  rep(0, length(tmax))
sum_vrn_fac <- sum_diff <- 0

#His vern function
vrn.fac = triangular(tavg,pbase,popt,pmax)
de.vrn = rep(0, length(tmax))
de.vrn[id$tmax > 30] = (tmax[tmax > 30] - 30) *   0.5

de.vrn[de.vrn > vrn.fac] = vrn.fac[de.vrn > vrn.fac]
ind = cumsum(vrn.fac - de.vrn) < 10
vrn.fac[ind] = vrn.fac[ind] - de.vrn[ind]
cum.vrn = cumsum(vrn.fac)
vrn.fac = cum.vrn/1
vrn.fac[vrn.fac > 1] = 1
vrn.fac[vrn.fac < 0] = 0

#My stan friendly version
for(i in 1:length(tavg)){
  
  vrn_fac <-  triangular(tavg[i],pbase,popt,pmax);
  de_vrn <- if_else(tmax[i] > 30,(tmax[i] - 30)/2,0)
  
  de_vrn <- if_else(de_vrn > vrn_fac,vrn_fac,de_vrn)
  sum_diff <- sum_diff + vrn_fac - de_vrn
  
  vrn_fac <- if_else(sum_diff<10,vrn_fac - de_vrn,vrn_fac)
  
  sum_vrn_fac <- sum_vrn_fac + vrn_fac
  vrn_fac <- max(min(sum_vrn_fac,1),0)
  
  final_vrn_fac[i] <- vrn_fac
  
}
################
#His Wheat Phenology function ignoring , from the "simplepheno" package
red_pheno <-function(tavg,doy,t.base,t.opt,t.max,tt.h,tt.hm){
    
    tt.daily = triangular(tavg,t.base,t.opt,t.max)*t.opt   
    dayl_fac <- ppsen_fun(doy,27.37177,70)
    tt.cum.adj = cumsum(tt.daily*dayl_fac)
    
    # Calculate days from sowing/germination to heading    
    dth = (1:length(tt.cum.adj))[tt.cum.adj>tt.h][1]
    
    tt.cum = cumsum(c(tt.daily[1:dth],tt.daily[(dth+1):length(tt.daily)]))
    
    # Calculate days from sowing/germination to maturity    
    dtm = (1:length(tt.cum))[tt.cum>(tt.h+tt.hm)][1]
    
    return(data.frame(dth,dtm))
    
}

stan_pheno <-function(tavg,t.base,t.opt,t.max,tt.h,tt.hm){
  
  tt.daily = stan_triangular(tavg,t.base,t.opt,t.max)*t.opt   
  tt.cum.adj = cumulative_sum(tt.daily)
  
  # Calculate days from sowing/germination to heading
  dth <- 1
  while(tt.cum.adj[dth]<tt.h){
    dth <- dth+1
  }
  
  dtm <- 1
  while(tt.cum.adj[dtm]<(tt.h+tt.hm)){
    dtm <- dtm+1
  }

  return(data.frame(dth,dtm))
  
}

library(simplepheno)
data(weather_temperate_2011)
data(p_triangular)
obstavg <- (weather_temperate_2011$tmin+weather_temperate_2011$tmax)/2
obsp <- list(base=0,opt=26,max=34)

doy_2011 <- as.POSIXlt(weather_temperate_2011$date)$yday+1
red_pheno(obstavg,doy_2011,0,26,34,tt.h=850,tt.hm=900)

stan_pheno(obstavg,0,25,30,tt.h=850,tt.hm=900)
