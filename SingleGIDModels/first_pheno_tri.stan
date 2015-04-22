functions{


  real stan_pheno(vector tavg, real tbase, real topt, real tmax, real tth){
    
    real ttcumadj;
    real z;
    int i;
    real daysto;
    real ttdaily;
    
    ttcumadj <- 0.0;
    daysto <- 1.0;
    i <- 1;
    while(ttcumadj<tth){

      z <- 1-step(tavg[i]-topt); 

      ttdaily <- (z*(tavg[i]-tbase)/(topt-tbase)+(1-z)*(tmax-tavg[i])/(tmax-topt));

      if(ttdaily<0){
        ttdaily <- 0;
      }

      ttdaily <- if_else(ttdaily>1,1,ttdaily*topt);

      ttcumadj <- ttcumadj+ttdaily;
      daysto <- daysto+1.0;
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
  real settth;

}

parameters {

  real theta1;
  real theta2;  
  real theta3;
  real<lower=0.5> sigma;

}

transformed parameters{

  real mu;
  mu <- stan_pheno(tavg, theta1, theta2, theta3, settth);

}

model {

  sigma ~ uniform(0.5,100);

  theta1 ~ uniform(tlower[1],tupper[1]);
  theta2 ~ uniform(tlower[2],tupper[2]);
  theta3 ~ uniform(tlower[3],tupper[3]);


  //print("dth=",mu);
  //print("obs=",dths);
  dths ~ normal(mu,sigma);

}