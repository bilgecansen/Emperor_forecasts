
data {
  int<lower=0> N;
  int<lower=0> K;
  int<lower=0> site_no[K];
  vector[N] X;
  real<lower=0> y[K];
  vector[N] y_mean;
  vector[K] y_sd;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> mu_sigma;
  real<lower=0> sigma_site;
  real<lower=0> sigma_time[N];
  real<lower=0> tau;
  real y_raw_time[K];
  real y_raw_site[N];
}

transformed parameters {
  real y_mu[N];
  real eps_site[N];
  real y_lat[K];
  
  for (i in 1:N) {
    eps_site[i] = y_raw_site[i]*sigma_site;
    y_mu[i] = X[i]*beta + alpha + eps_site[i];
  }
  
  for (i in 1:K) {
    y_lat[i] = y_mu[site_no[i]] + y_raw_time[i]*sigma_time[site_no[i]];
  }
}

model {
  
  y_raw_site ~ normal(0, 1);
  y_raw_time ~ normal(0, 1);
  mu_sigma ~ normal(0, 5);
  tau ~ normal(0, 5);
  alpha ~ normal(0, 20);
  beta ~ normal(0, 10);
  sigma_site ~ normal(0,5);
  
  sigma_time ~ normal(mu_sigma, tau);
  
  for (i in 1:K){
    y[i] ~ normal(y_lat[i], y_sd[i]);
  } 
}

generated quantities {
  real Rsq;
  vector[N] res;
  vector[N] error;
  vector[N] mu;
  
  mu = X*beta + alpha;
  
  for (i in 1:N) {
    res[i] = (mu[i]- mean(y_mean))^2;
    error[i] = (y_mean[i] - mu[i])^2;
  }
  
  Rsq = (sum(res)/(N-1))/((sum(res)/(N-1)) + (sum(error)/(N-1)));
  
  int<lower = 0, upper = 1> mean_gt;
  int<lower = 0, upper = 1> sd_gt;
  vector[K] y_rep;
  
  for (i in 1:K) {
    y_rep[i] = normal_rng(y_lat[i], y_sd[i]);
  }
  
  mean_gt = mean(y_rep) > mean(y);
  sd_gt = sd(y_rep) > sd(y);
  
}
