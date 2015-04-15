functions{


  real stan_pheno(vector tavg, real tbase, real topt, real tmax, real tth){
    
    real ttcumadj;
    real z;
    int i;
    real daysto;
    real ttdaily;
    real a;
    real denom;

    ttcumadj <- 0.0;
    daysto <- 0.0;
    i <- 1;
    
    //Constants used in fitting the wang function
    a  <-  log(2) / log((tmax - tbase) / (topt - tbase));
    denom <- pow(topt - tbase,2*a);

    while(ttcumadj<tth){

      daysto <- daysto+1.0;

      ttdaily <- (2*pow(tavg[i] - tbase,a)*pow(topt - tbase,a) - pow(tavg[i] - tbase,2*a)) / denom;
      ttdaily <- if_else(tavg[i]<tbase||tavg[i]>tmax,0,ttdaily*topt);

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
  vector[3] tlower;
  vector[3] tupper;
  real tthLow;
  real tthHigh;

}

parameters {

  real theta1;
  real theta2;  
  real theta3;
  real<lower=0> sigma;
  real tthpar;

}

transformed parameters{

  real mu;
  mu <- stan_pheno(tavg, theta1, theta2, theta3, tthpar);

}

model {

  sigma ~ uniform(0,10);

  theta1 ~ uniform(tlower[1],tupper[1]);
  theta2 ~ uniform(tlower[2],tupper[2]);
  theta3 ~ uniform(tlower[3],tupper[3]);

  tthpar ~ uniform(tthLow,tthHigh);
  dths ~ normal(mu,sigma);

}