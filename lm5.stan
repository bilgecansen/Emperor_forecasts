
data {
  int<lower=0> N;
  int<lower=0> K;
  //int<lower=0> T;
  int<lower=0> site_no[K];
  int<lower=0> year_no[K];
  vector[N] X;
  real<lower=0> y[K];
  vector[K] y_sd;
}

parameters {
  real alpha;
  real beta;
  //real<lower=0> mu_sigma;
  real<lower=0> sigma_site;
  real<lower=0> sigma_time[N];
  //real<lower=0> tau;
  //real sigma_raw_time[N];
  real y_raw_time[K];
  real y_raw_site[N];
}

transformed parameters {
  real y_mu[N];
  //real<lower=0> sigma_time[N];
  real eps_site[N];
  real y_lat[K];
  
  for (i in 1:N) {
    //sigma_time[i] = mu_sigma + sigma_raw_time[i]*tau;
    eps_site[i] = y_raw_site[i]*sigma_site;
    y_mu[i] = X[i]*beta + alpha + eps_site[i];
  }
  
  for (i in 1:K) {
    y_lat[i] = y_mu[site_no[i]] + y_raw_time[i]*sigma_time[site_no[i]];
  }
}

model {
  
  //sigma_raw_time ~ normal(0, 1);
  y_raw_site ~ normal(0, 1);
  y_raw_time ~ normal(0, 1);
  //mu_sigma ~ normal(0, 1);
  //tau ~ normal(0, 1);
  sigma_time ~ normal(0,0.5);
  alpha ~ normal(9, 2);
  beta ~ normal(0, 1);
  
  for (i in 1:K){
    y[i] ~ normal(y_lat[i], y_sd[i]);
  } 
}
