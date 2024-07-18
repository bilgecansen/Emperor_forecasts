
data {
  int<lower=0> N;
  int<lower=0> T;
  int<lower=0> K;
  int<lower=0> site_no[K];
  int<lower=0> year_no[K];
  vector[N] X;
  vector[K] y;
  vector[K] y_sd;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> mu;
  real<lower=0> epsilon;
  real<lower=0> sigma;
  real<lower=0> sigma_time;
  //real y_raw[N];
  real y_raw_time[K];
}

transformed parameters {
  real y_lat[K];
  //real eps[N];
  real eps_time[K];
  
  //for (i in 1:N) {
    //eps[i] = y_raw[i]*sigma;
    //for (t in 1:T) {
      //eps_time[i,t] = y_raw_time[i,t]*sigma_time;
      //y_lat[i,t] = X[i,t]*beta + alpha + eps_time[i,t];
    //}
  //}
  
  for (i in 1:K) {
    eps_time[i] = y_raw_time[i]*sigma_time;
    y_lat[i] = X[site_no[i]]*beta + alpha + eps_time[i];
  }
}

model {
  sigma ~ normal(0, 1);
  mu ~ normal(0, 1);
  epsilon ~ normal(0, 1);
  sigma_time ~ normal(0, 1);
  y_raw_time ~ normal(0, 1);
  //for (i in 1:N) y_raw_time[i,] ~ normal(0, 1);
  
  for (i in 1:K){
    y[i] ~ normal(y_lat[i], y_sd[i]);
  } 
}
