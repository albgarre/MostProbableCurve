

data {
  int<lower=0> n;
  int count[n];
  vector[n] t;
  int dil[n];
  vector[n] vol;
  vector[n] temperature;
  real Tref;
  
}


parameters {
  real logN0;
  real logDref;
  real z;
}

transformed parameters {

  real<lower=0> Dref;    
  real<lower=0> N0;
  Dref = 10^logDref;
  N0 = 10^logN0;
}

model {
  
  // Likelihood
  vector[n] lambda;
  vector[n] Nmedia;
  vector[n] logD;
  vector[n] D;
  for (i in 1:n)
    logD[i] = logDref - (temperature[i] - Tref)/z;
  for (i in 1:n)
    D[i] = 10^logD[i];
  for (i in 1:n)
    Nmedia[i] = N0*10^(-t[i]/D[i]);
  for (i in 1:n)
    lambda[i] = Nmedia[i]*vol[i]/(10^dil[i]);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
  count ~ poisson(lambda);
  
  // Priors
  
  logD ~ normal(0, 2);
  logN0 ~ normal(8, 100);

}










