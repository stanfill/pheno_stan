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

################
#His Wheat Phenology function ignoring , from the "simplepheno" package
red_pheno <-function(tavg,t.base,t.opt,t.max,tt.h,tt.hm){
    
    tt.daily = triangular(tavg,t.base,t.opt,t.max)*t.opt   
    tt.cum.adj = cumsum(tt.daily)
    
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

red_pheno(obstavg,0,26,34,tt.h=850,tt.hm=900)
stan_pheno(obstavg,0,25,30,tt.h=850,tt.hm=900)
