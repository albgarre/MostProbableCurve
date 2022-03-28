
data {
  int<lower=0> n;
  vector[n] logN;
  vector[n] t;
}


parameters {
  real logN0;
  real logD;
  real logp;
  real<lower=0> sigma;
}

transformed parameters {

  real<lower=0> D;    
  real<lower=0> N0;
  real<lower=0> p;
  D = 10^logD;
  N0 = 10^logN0;
  p = 10^logp;
}

model {
  
  // Likelihood
  vector[n] E;
  for (i in 1:n)
    E[i] = logN0 - (t[i]/D)^p;
  
  for (i in 1:n)
    logN[i] ~ normal(E[i], sigma);
  
  // Priors
  
  logD ~ normal(0, 2);
  logp ~ normal(0, 1);
  logN0 ~ normal(8, 100);
  sigma ~ exponential(1);

}


