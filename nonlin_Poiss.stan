
data {
  int<lower=0> n;
  int count[n];
  vector[n] t;
}


parameters {
  real logN0;
  real logdelta;
  real logp;
}

transformed parameters {

  real<lower=0> p;
  real<lower=0> delta;    
  real<lower=0> N0;
  p = 10^logp;
  delta = 10^logdelta;
  N0 = 10^logN0;
}

model {
  
  // Likelihood
  vector[n] lambda;
  for (i in 1:n)
    lambda[i] = N0*10^(-(t[i]/delta)^p);

  count ~ poisson(lambda);
  
  // Priors
  
  logp ~ normal(0, 1);
  logdelta ~ normal(0, 2);
  logN0 ~ normal(8, 100);

}







