functions{


  vector stan_pheno(row_vector tavg, real tbase, real topt, real tmax, real tth, real tthm){
    
    real ttcumadj;
    real z;
    int i;
    vector[2] daysto;
    real ttdaily;
    real ttm;
    
    ttm <- tth+tthm;

    ttcumadj <- 0.0;
    daysto[1] <- 0.0;
    daysto[2] <- 0.0;

    i <- 1;
    while(ttcumadj<ttm){

      daysto[1] <- if_else(ttcumadj<tth,daysto[1]+1.0,daysto[1]);
      daysto[2] <- daysto[2]+1.0;

      //This is the triangular phenology model
      z <- 1-step(tavg[i]-topt); 
      ttdaily <- (z*(tavg[i]-tbase)/(topt-tbase)+(1-z)*(tmax-tavg[i])/(tmax-topt));

      if(ttdaily<0){
        ttdaily <- 0;
      }

      ttdaily <- if_else(ttdaily>1,1,ttdaily*topt);

      //End the triangular model      

      ttcumadj <- ttcumadj+ttdaily;

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

  matrix[nyears,ndays] tavg;
  real obs_dth[nyears,ngid,nobs];
  real obs_dtm[nyears,ngid,nobs];
  
  vector[3] tlower;
  vector[3] tupper;
  real tthLow;
  real tthHigh;
  real tthmLow;
  real tthmHigh;

}

parameters {
  
  real tmin; //Phenology model parameters
  real topt;  
  real tmax;

  real<lower=0> sigma_dth;      //Residual variance
  real<lower=0> sigma_dtm;      //Residual variance


  real tth_g[ngid];        //Genome specific tth value
  real tthm_g[ngid];        //Genome specific tth value


  real mu_tth;              //Mean of the tthpars
  real<lower=0> sig_tth;    //sd of tthpars
  real mu_tthm;              //Mean of the ttmpars
  real<lower=0> sig_tthm;    //sd of ttmpars
}

transformed parameters{

  real dth_hat_g[nyears,ngid];
  real dtm_hat_g[nyears,ngid];
  vector[2] mulk;

  for(l in 1:nyears){
    for(k in 1:ngid){ 
      mulk <- stan_pheno(tavg[l], tmin, topt, tmax, tth_g[k],tthm_g[k]);
      dth_hat_g[l,k] <- mulk[1];
      dtm_hat_g[l,k] <- mulk[2];
    }
  }

}

model {

  //Hierarchy structure for tth parameter
  mu_tth ~ uniform(tthLow,tthHigh);
  sig_tth ~ uniform(0,5);
  tth_g ~ normal(mu_tth,sig_tth);

  //Hierarchy structure for tthm parameter
  mu_tthm ~ uniform(tthmLow,tthmHigh);
  sig_tthm ~ uniform(0,5);
  tthm_g ~ normal(mu_tthm,sig_tthm);

  sigma_dth ~ uniform(0,10);
  sigma_dtm ~ uniform(0,10);

  tmin ~ uniform(tlower[1],tupper[1]);
  topt ~ uniform(tlower[2],tupper[2]);
  tmax ~ uniform(tlower[3],tupper[3]);
  
  for(l in 1:nyears){
    for(n in 1:ngid){
      obs_dth[l,n] ~ normal(dth_hat_g[l,n],sigma_dth);
      obs_dtm[l,n] ~ normal(dtm_hat_g[l,n],sigma_dtm);
    }
  }

}