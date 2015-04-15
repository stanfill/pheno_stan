functions{


  real stan_pheno(row_vector tavg, real tbase, real topt, real tmax, real tth){
    
    real ttcumadj;
    real z;
    int i;
    real daysto;
    real ttdaily;
    
    ttcumadj <- 0.0;
    daysto <- 1.0;
    i <- 1;
    while(ttcumadj<tth){

      //This is the triangular phenology model
      z <- 1-step(tavg[i]-topt); 
      ttdaily <- (z*(tavg[i]-tbase)/(topt-tbase)+(1-z)*(tmax-tavg[i])/(tmax-topt));

      if(ttdaily<0){
        ttdaily <- 0;
      }

      ttdaily <- if_else(ttdaily>1,1,ttdaily*topt);
      //End the triangular model      

      ttcumadj <- ttcumadj+ttdaily;
      daysto <- daysto+1.0;
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
  real dths[nyears,ngid,nobs];
  
  vector[3] tlower;
  vector[3] tupper;
  real tthLow;
  real tthHigh;

}

parameters {
  
  real theta1; //Phenology model parameters
  real theta2;  
  real theta3;

  real<lower=0> sigma;      //Residual variance

  real tthpar[ngid];        //Genome specific tth value

  real mu_tth;              //Mean of the tthpars
  real<lower=0> sig_tth;    //sd of tthpars
}

transformed parameters{

  real mu[nyears,ngid];

  for(l in 1:nyears){
    for(k in 1:ngid){ 
      mu[l,k] <- stan_pheno(tavg[l], theta1, theta2, theta3, tthpar[k]);
    }
  }

}

model {

  mu_tth ~ uniform(tthLow,tthHigh);
  sig_tth ~ uniform(0,5);

  tthpar ~ normal(mu_tth,sig_tth);

  sigma ~ uniform(0,10);

  theta1 ~ uniform(tlower[1],tupper[1]);
  theta2 ~ uniform(tlower[2],tupper[2]);
  theta3 ~ uniform(tlower[3],tupper[3]);
  
  for(l in 1:nyears)
    for(n in 1:ngid)
      dths[l,n] ~ normal(mu[l,n],sigma);

}