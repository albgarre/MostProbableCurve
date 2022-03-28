
data {
  int<lower=0> n;
  int count[n];
  vector[n] t;
}


parameters {
  real logN0;
  real logD;
}

transformed parameters {

  real<lower=0> D;    
  real<lower=0> N0;
  D = 10^logD;
  N0 = 10^logN0;
}

model {
  
  // Likelihood
  vector[n] lambda;
  for (i in 1:n)
    lambda[i] = N0*10^(-t[i]/D);

  count ~ poisson(lambda);
  
  // Priors
  
  logD ~ normal(0, 2);
  logN0 ~ normal(8, 1);

}







