

data {
  int<lower=0> n;
  int<lower=1> nbio;
  int count[n];
  vector[n] t;
  int dil[n];
  vector[n] vol;
  int bio[n];
  
}


parameters {
  // log N0
  real logN0[nbio];
  real m_logN0;
  real s_logN0;
  
  // log D
  real logD[nbio];
  real m_logD;
  real s_logD;
}

transformed parameters {

  real<lower=0> D[nbio];    
  real<lower=0> N0[nbio];
  
  for (i in 1:nbio)
    N0[i] = 10^logN0[i];
    
  for (i in 1:nbio)
    D[i] = 10^logD[i];
}

model {

  // Likelihood
  vector[n] lambda;
  vector[n] Nmedia;
  
  for (i in 1:n)
    Nmedia[i] = N0[bio[i]]*10^(-t[i]/D[bio[i]]);
  for (i in 1:n)
    lambda[i] = Nmedia[i]*vol[i]/(10^dil[i]);
  
  // Likelihood 
  
  count ~ poisson(lambda);
  
  // Variability level
  
  logN0 ~ normal(m_logN0, s_logN0);
  logD ~ normal(m_logD, s_logD);
  
  // Priors
  
  // logD ~ normal(0, 2);
  // logN0 ~ normal(8, 100);
  
  m_logN0 ~ normal(8, 100);
  s_logN0 ~ exponential(2);
  
  m_logD ~ normal(0, 2);
  s_logD ~ exponential(1);

}










