functions{

  real calc_dayl_fac(real lat, real doy, real ppsen){
  
    real s1;
    real c1;
    real dec;
    real dlv;
    real dayl;
    real dayl_fac;

    s1  <-  sin(lat * 0.01745);
    c1  <-  cos(lat * 0.01745);
    dec  <-  0.4093 * sin(0.0172 * (doy - 82.2));
    dlv  <-  (((-s1) * sin(dec) - 0.1047)/(c1 * cos(dec)));
  
    if(dlv < (-0.87)){
      dlv <- (-0.87);
    }

    dayl  <-  7.639 * acos(dlv);
    dayl_fac  <-  1 - (ppsen/10000) * pow(20 - dayl,2);
  
    return dayl_fac;

  }

  real wang_pheno(real temp, real tbase, real topt, real tmax){

    real a;
    real denom;
    real ttdaily;
  
    //Constants used in fitting the Wang function
    a  <-  (log2()) / log((tmax - tbase) / (topt - tbase));
    denom <- pow(topt - tbase,2*a);
      
    //This is the Wang (curved) phenology model
    ttdaily <- (2*pow(temp - tbase,a)*pow(topt - tbase,a) - pow(temp - tbase,2*a)) / denom;
  
    if(temp<tbase||temp>tmax)
      ttdaily <- 0;
    else
      ttdaily <- ttdaily*topt;

    return ttdaily;
      
  }

  real triangle_pheno(real temp, real tbase, real topt, real tmax){

    real z;
    real ttdaily;
  
    z <- 1-step(temp-topt); 
    ttdaily <- (z*(temp-tbase)/(topt-tbase)+(1-z)*(tmax-temp)/(tmax-topt));

    if(ttdaily<0)
      ttdaily <- 0;

    if(ttdaily>1)
      ttdaily <- 1;
    
    ttdaily <- ttdaily*topt;

    return ttdaily;
      
  }


  vector stan_pheno(row_vector tavg, row_vector doy, row_vector obs_tmax, real tbase, real topt, real tmax, real tth, real tthm, real ppsen){
    
    int n_obs;
    real ttcumadj;
    real z;
    int i;
    vector[2] daysto;
    real ttdaily;
    real ttm;
  
    real dayl_fac;

    n_obs <- num_elements(tavg); //To protect against going too far in tavg vector

    ttm <- tth+tthm;
    ttcumadj <- 0.0;
    daysto[1] <- 0.0;
    daysto[2] <- 0.0;

    i <- 1;

    while(ttcumadj<ttm && i<n_obs){

      daysto[2] <- daysto[2]+1.0;

      //Use the Triangle phenology model to calculate day i thermal time
      ttdaily <- triangle_pheno(tavg[i],tbase,topt,tmax);


      if(ttcumadj<tth){
  
        //Calculate day lengt factor at lat=27.37177, ppsen variable
        dayl_fac <- calc_dayl_fac(27.37177, doy[i], ppsen);

        //Until heading is reached, adjust the daily tt for dayl and vern_fac
        daysto[1] <- daysto[1]+1.0;
        ttcumadj <- ttcumadj+ttdaily*dayl_fac;
      }
      else{
        //After heading is reached, don't make any adjustments
        ttcumadj <- ttcumadj+ttdaily;
      }

      i <- i+1;
    }

    return daysto;
  
  }

}

data {
  
  int<lower=0> ndays;
  int<lower=0> nobs;
  int<lower=0> ngid;
  int<lower=0> nyears;

  matrix[nyears,ndays] obs_tavg;
  matrix[nyears,ndays] doy;
  matrix[nyears,ndays] obs_tmax;

  real obs_dth[nyears,ngid,nobs];
  real obs_dtm[nyears,ngid,nobs];
  
  real tthLow;
  real tthHigh;
  real tthmLow;
  real tthmHigh;

  //real ppsen;
  real topt;

  real tmin;
  real tmax;

}

parameters {
  
  //real<lower=20,upper=30> topt;
  real ppsen[ngid];
  real<lower=0> sigma_dth;      //Residual variance
  real<lower=0> sigma_dtm;      //Residual variance


  real<lower=600,upper=1400> tth_g[ngid];        //Genome specific tth value
  real<lower=600,upper=1400> tthm_g[ngid];        //Genome specific tth value

}


model {

  //I don't actually care about the current estimate of dth or dtm so make them 
  //local variables that don't exist outside this code chunck
  vector[2] mulk;

  tth_g ~ normal(tthLow,tthHigh);
  tthm_g ~ normal(tthmLow,tthmLow);              //Each tthm is assumed normal

  sigma_dth ~ uniform(0,40);
  sigma_dtm ~ uniform(0,70);

  
  ppsen ~ uniform(30, 95);
  //topt ~ uniform(20,30);

  for(l in 1:nyears){
    for(n in 1:ngid){

      mulk <- stan_pheno(obs_tavg[l], doy[l], obs_tmax[l], tmin, topt, tmax, tth_g[n],tthm_g[n], ppsen[n]);

      obs_dth[l,n] ~ normal(mulk[1],sigma_dth);
      obs_dtm[l,n] ~ normal(mulk[2],sigma_dtm);

    }
  }

}
