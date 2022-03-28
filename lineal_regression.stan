
data {
  int<lower=0> n;
  vector[n] logN;
  vector[n] t;
}


parameters {
  real logN0;
  real logD;
  real<lower=0> sigma;
}

transformed parameters {

  real<lower=0> D;    
  real<lower=0> N0;
  D = 10^logD;
  N0 = 10^logN0;
}

model {
  
  // Likelihood
  vector[n] E;
  for (i in 1:n)
    E[i] = logN0 - t[i]/D;
  
  for (i in 1:n)
    logN[i] ~ normal(E[i], sigma);
  
  // Priors
  
  logD ~ normal(0, 2);
  logN0 ~ normal(8, 2);
  sigma ~ exponential(1);

}


