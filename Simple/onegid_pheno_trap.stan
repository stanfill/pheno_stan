functions{


  real stan_pheno(vector tavg, real tbase, real topt1, real topt2, real tmax, real tth){
    
    real ttcumadj;
    real z;
    int i;
    real daysto;
    real ttdaily;
    real toptmid;
  
    toptmid <- (topt1+topt1)/2;
    ttcumadj <- 0.0;
    daysto <- 0.0;
    i <- 1;

    while(ttcumadj<tth){

      daysto <- daysto+1.0;

      ttdaily <- if_else(tavg[i]<=topt1,(tavg[i] - tbase) / (topt1 - tbase),0);
      ttdaily <- if_else(tavg[i]>topt2,(tmax - tavg[i]) / (tmax-topt2),0);
      if(ttdaily<0){
        ttdaily <- 0;
      }
      ttdaily <- if_else(ttdaily>1,1,ttdaily*toptmid);


      ttcumadj <- ttcumadj+ttdaily;
      i <- i+1;
    }


    return daysto;
  
  }

}

data {
  
  int<lower=1> ndays;
  int<lower=1> nobs;
  vector[ndays] tavg;
  vector[nobs] dths;
  vector[4] tlower;
  vector[4] tupper;
  real tthLow;
  real tthHigh;

}

parameters {

  real theta1;
  real theta2;  
  real theta3;
  real theta4;
  real<lower=0> sigma;
  real tthpar;

}

transformed parameters{

  real mu;
  mu <- stan_pheno(tavg, theta1, theta2, theta3, theta4, tthpar);

}

model {

  sigma ~ uniform(0,10);

  theta1 ~ uniform(tlower[1],tupper[1]);
  theta2 ~ uniform(tlower[2],tupper[2]);
  theta3 ~ uniform(tlower[3],tupper[3]);
  theta4 ~ uniform(tlower[4],tupper[4]);

  tthpar ~ uniform(tthLow,tthHigh);

  //print("dth=",mu);
  //print("obs=",dths);
  dths ~ normal(mu,sigma);

}