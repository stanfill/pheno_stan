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
  
  real theta1; //Phenology model parameters
  real theta2;  
  real theta3;

  real<lower=0> tth_resid;      //Residual variance
  real<lower=0> tthm_resid;      //Residual variance


  real tthpar[ngid];        //Genome specific tth value
  real tthmpar[ngid];        //Genome specific tth value


  real mu_tth;              //Mean of the tthpars
  real<lower=0> sig_tth;    //sd of tthpars
  real mu_tthm;              //Mean of the ttmpars
  real<lower=0> sig_tthm;    //sd of ttmpars
}

transformed parameters{

  real mutth_g[nyears,ngid];
  real muttm_g[nyears,ngid];
  vector[2] mulk;

  for(l in 1:nyears){
    for(k in 1:ngid){ 
      mulk <- stan_pheno(tavg[l], theta1, theta2, theta3, tthpar[k],tthmpar[k]);
      mutth_g[l,k] <- mulk[1];
      muttm_g[l,k] <- mulk[2];
    }
  }

}

model {

  //Hierarchy structure for tth parameter
  mu_tth ~ uniform(tthLow,tthHigh);
  sig_tth ~ uniform(0,5);
  tthpar ~ normal(mu_tth,sig_tth);

  //Hierarchy structure for tthm parameter
  mu_tthm ~ uniform(tthmLow,tthmHigh);
  sig_tthm ~ uniform(0,5);
  tthmpar ~ normal(mu_tthm,sig_tthm);

  tth_resid ~ uniform(0,10);
  tthm_resid ~ uniform(0,10);

  theta1 ~ uniform(tlower[1],tupper[1]);
  theta2 ~ uniform(tlower[2],tupper[2]);
  theta3 ~ uniform(tlower[3],tupper[3]);
  
  for(l in 1:nyears){
    for(n in 1:ngid){
      obs_dth[l,n] ~ normal(mutth_g[l,n],tth_resid);
      obs_dtm[l,n] ~ normal(muttm_g[l,n],tthm_resid);
    }
  }

}