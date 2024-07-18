
data {
  int<lower=0> N;
  int<lower=0> L;
  vector[N] X;
  vector[L] X_pred;
  vector[N] y;
  vector[N] y_sd;
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
  real y_raw[N];
}

transformed parameters {
  real y_lat[N];
  real eps[N];
  
  for (i in 1:N) {
    eps[i] = y_raw[i]*sigma;
    y_lat[i] = X[i]*beta + alpha + eps[i];
  }
  
}

model {
  sigma ~ normal(0, 5);
  y_raw ~ normal(0, 1);
  y ~ normal(y_lat, y_sd);
}

generated quantities {
  real Rsq;
  vector[N] res;
  vector[N] error;
  vector[L] mu_pred;
  vector[N] mu;
  
  mu = X*beta + alpha;
  mu_pred = X_pred*beta + alpha;
  
  for (i in 1:N) {
    res[i] = (mu[i]- mean(y))^2;
    error[i] = (y[i] - mu[i])^2;
  }
  
  Rsq = (sum(res)/(N-1))/((sum(res)/(N-1)) + (sum(error)/(N-1)));
}
