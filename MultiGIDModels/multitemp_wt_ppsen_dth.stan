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
    else
      ttdaily <- ttdaily*topt;

    return ttdaily;
      
  }

  
  vector calc_vern_sens(real tavg, real obs_tmax, real pbase, real popt, real pmax, real sum_diff){

    //Calculate the vernalisation factor for a given day.  Vernalisation is cumulative
    //so it needs to be added as it goes.  The first element of the returned vector (vrn_fac_etal)
    //is the vernalisation factor for day i.  The second element is the cumulative sum
    //of vrn_fac-de_vrn

    real vrn_fac;
    real de_vrn;
    real sum_vrn_fac;
    vector[2] vrn_fac_etal;
    real diff_cumsum;

    de_vrn <- 0.0;    

    vrn_fac <-  triangle_pheno(tavg,pbase,popt,pmax);
    vrn_fac <- vrn_fac/popt;  

    if(obs_tmax > 30)
      de_vrn <- (obs_tmax - 30)/2;
  
    if(de_vrn > vrn_fac)
      de_vrn <- vrn_fac;
    
    diff_cumsum <- sum_diff + vrn_fac - de_vrn;
  
    if(diff_cumsum<10)
      vrn_fac <- vrn_fac - de_vrn;

    vrn_fac_etal[1] <- vrn_fac;
    vrn_fac_etal[2] <- diff_cumsum;

    return vrn_fac_etal;  
  }

  real stan_pheno(row_vector tavg, row_vector doy, row_vector obs_tmax, real tbase, real topt, real tmax, real tth, real ppsen){
    
    int n_obs;
    real ttcumadj;
    real z;
    int i;
    real daystoh;
    vector[2] ver_fac_res;
    real ttdaily;
  
    real dayl_fac;
    real vern_fac;
    real vf_i;

    n_obs <- num_elements(tavg); //To protect against going too far in tavg vector

    vern_fac <- 0.0;
    ver_fac_res[2] <- 0.0;
    ttcumadj <- 0.0;
    daystoh <- 0.0;


    i <- 1;

    while(ttcumadj<tth && i<n_obs){

      //Use the Triangle phenology model to calculate day i thermal time
      ttdaily <- triangle_pheno(tavg[i],tbase,topt,tmax);

      //Calculate day lengt factor at lat=27.37177, ppsen=30
      dayl_fac <- calc_dayl_fac(27.37177, doy[i], ppsen);

      //Calculate vernalisation factor with pbase=-5, popt=7 and pmax=15
      //ver_fac_res <- calc_vern_sens(tavg[i], obs_tmax[i], (-5.0), 7.0, 15.0,ver_fac_res[2]);
      ver_fac_res <- calc_vern_sens(tavg[i], obs_tmax[i], 0.0, 26.0, 34.0,ver_fac_res[2]);
      vf_i <- ver_fac_res[1];
      vern_fac <- vern_fac+vf_i;

      if(vern_fac>1)
        vern_fac <- 1;
      
      if(vern_fac<0)
        vern_fac <- 0;


      //Until heading is reached, adjust the daily tt for dayl and vern_fac
      daystoh <- daystoh +1.0;
      ttcumadj <- ttcumadj+ttdaily*dayl_fac*vern_fac;

      i <- i+1;
    }

    return daystoh;
  
  }

}

data {
  
  int<lower=0> ndays;
  int<lower=0> nobs;
  int<lower=0> nyears;

  matrix[nyears,ndays] obs_tavg;
  matrix[nyears,ndays] doy;
  matrix[nyears,ndays] obs_tmax;

  matrix[nyears,nobs] obs_dth;
  
  vector[3] tlower;
  vector[3] tupper;
  real tthLow;
  real tthHigh;


}

parameters {
  
  real<lower=(-5),upper=5> tmin; //Phenology model parameters
  real<lower=20,upper=30> topt;  
  real<lower=30,upper=50> tmax;
  real ppsen;


  real<lower=0> sigma_dth;      //Residual variance


  real<lower=700,upper=1400> tth_g[nyears];        //Genome specific tth value


  //real<lower=700, upper=1400> mu_tth;              //Mean of the tthpars
  //real<lower=0> sig_tth;    //sd of tthpars

}


model {

  //I don't actually care about the current estimate of dth or dtm so make them 
  //local variables that don't exist outside this code chunck
  real mulk;

  //Hierarchy structure for tth parameter
  //mu_tth ~ normal(tthLow,tthHigh);
  //sig_tth ~ uniform(0,5);
  tth_g ~ normal(tthLow,tthHigh);

  sigma_dth ~ uniform(0,40);

  tmin ~ normal(tlower[1],tupper[1]);
  topt ~ normal(tlower[2],tupper[2]);
  tmax ~ normal(tlower[3],tupper[3]);
  
  ppsen ~ normal(30, 5);


  for(l in 1:nyears){

      mulk <- stan_pheno(obs_tavg[l], doy[l], obs_tmax[l], tmin, topt, tmax, tth_g[l], ppsen);
      //print("mulk = ", mulk); //Compare to results of pheno()
      obs_dth[l] ~ normal(mulk,sigma_dth);

    
  }
  //print("nyears = ",nyears);
}